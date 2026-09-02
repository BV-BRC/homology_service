# Rebuilding the BLAST databases

The homology service searches a *generation*: a dated directory of BLAST
databases plus a sqlite catalog that maps taxon ids onto them. `deploy.cfg`
names the live one in two places. Rebuilding means constructing a new
generation alongside the old one and repointing at the end, never editing one
in place.

A generation looks like this:

    data.<YYYY-MMDD>/
        by-genus-bacterial/     one database per genus     ~11,970 databases
        by-genus-viral/         one per *family*, despite the directory name    ~195
        ref/                    curated reference sets     6 databases
        db.sqlite               the catalog
        logs/

Each database is really three, one per `<dbtype> <ftype>` pair, and each is
accompanied by `.taxids`, `.taxlist` and `.glist` sidecars:

| built as | file | contents |
|---|---|---|
| `aa features` | `<name>.features.faa` | protein features |
| `dna features` | `<name>.features.fna` | DNA features |
| `dna contigs` | `<name>.contigs.fna` | contigs |

## Before you start

**Disk.** A full generation is about **1.8 TB**, essentially all of it
bacterial: 703 GB of `features.faa`, 641 GB of `features.fna`, 415 GB of
`contigs.fna`, plus ~50 GB for viral and `ref`. Confirm the target has that
much *before* the first invocation — the failure mode is not a clean stop, see
"When the disk fills" below.

Build on a local filesystem rather than on `/vol/blastdb`. The reason is
contention, not speed: `/vol/blastdb` is a shared netapp serving live BLAST
queries, and 1.8 TB of writes degrades production and is degraded by it. The
build's own I/O is large-block sequential, so local spinning disk versus SSD is
a second-order difference — the one place SSD genuinely helps is that 12 build
threads times `--parallel 2` means up to 24 interleaved append streams. The
cost of building locally is a ~1.8 TB copy to `/vol/blastdb` at publish time.

**Perl modules.** `p3x-create-blast-db` loads `PerlIO::via::Blockwise` (which
needs `IO::AIO`) and `BerkeleyDB`. Neither is loaded by the driver, so a
missing one lets `p3x-build-pathogen-blastdbs` start normally and kills every
worker at compile time. The symptom is opaque:

    Failure running p3x-create-blast-db ...: No such file or directory 512

`$?` of 512 is exit 2, and Perl's `die` exits with `$!` when it is nonzero, so
`Can't locate Foo.pm in @INC` (ENOENT) arrives as exit 2 with the parent's
stale `$!` printed instead of the child's actual message. The real error is in
the child's log, and if the child died at compile time there is no log at all —
no log files anywhere is itself the diagnostic. One-second check:

    perl -c scripts/p3x-create-blast-db.pl

**BLAST binaries.** `Bio::P3::HomologySearch::HomologySearch` initializes PATH
on import, looking under `$KB_TOP/services/homology_service/bin` and
`$KB_TOP/modules/homology_service/blast.bin`. The JSON output format this code
relies on is version-specific; the tree pins BLAST+ 2.13.

## Stage 0 — a fresh dated directory

    export NEW=/blast-build/blast-$(date +%Y-%m%d)
    mkdir -p $NEW/{by-genus-bacterial,by-genus-viral,ref,logs}

**This is a correctness requirement, not tidiness.** `p3x-create-blast-db` skips
any database that already looks built:

```perl
if (-s $taxids && @blast_files && !$opt->overwrite)
{
    print STDERR "Blast data in $dbfile already exists, skipping build\n";
    exit 0;
}
```

`p3x-build-pathogen-blastdbs` grew an `--overwrite` (`-f`) passthrough in
2026-09, so a rerun *can* now be forced -- but for a full generation the fresh
directory is still the rule, because forcing an overwrite across a populated
tree rebuilds all 1.8 TB rather than the part you meant. Without the flag a
rerun into a populated directory cannot rebuild anything, and that is how
`data.2022-0916` shipped damaged: a rerun on 09-20 was refused for 5,110
databases whose logs contain nothing but "already exists, skipping build",
leaving the broken 09-18 copies in place. Confirm as you go:

    grep -rl 'already exists, skipping' $NEW/logs/ | wc -l    # must stay 0

For a *targeted repair* of an existing generation the opposite holds -- skips
are the point, and that count is meaningless. See "Recovering a partial build".

## Stage 1 — bacterial, by genus

Three passes over every genus, one per database type:

```bash
for spec in "aa features" "dna features" "dna contigs"; do
    set -- $spec
    p3x-build-pathogen-blastdbs --all-genera \
        --n-build-threads 12 --n-db-threads 2 \
        --log-dir $NEW/logs/bacterial.$2.$1 \
        $1 $2 $NEW/by-genus-bacterial
done
```

`--all-genera` discards the script's hardcoded 27-genus pathogen list and
queries the taxonomy for every genus instead. Work is distributed by
`LPTScheduler`, which bin-packs genera longest-first, so the big ones
(Bacillus, Streptococcus, Pseudomonas, Streptomyces) start immediately and the
tail is thousands of tiny genera.

`--n-build-threads` is the number of concurrent genera;
`--n-db-threads` becomes `--parallel` on each `p3x-create-blast-db`, so peak
concurrency is the product. `--build-tempdir` (default `/dev/shm`) is used only
for genera under 50,000 in size; larger ones get `TMPDIR=/disks/tmp`, which the
driver sets unconditionally at startup.

Per-genus stdout and stderr land in `$NEW/logs/bacterial.<ftype>.<dbtype>/log.<Genus>.{out,err}`.
Every genus produces a non-empty `.err` in normal operation — it carries
"Skipping duplicate id" and "Bad sequence" warnings — so a non-empty `.err` is
not a failure signal by itself.

## Stage 2 — viral, by family

```bash
for spec in "aa features" "dna features" "dna contigs"; do
    set -- $spec
    p3x-build-pathogen-blastdbs --viral \
        --n-build-threads 12 --n-db-threads 2 \
        --log-dir $NEW/logs/viral.$2.$1 \
        $1 $2 $NEW/by-genus-viral
done
```

`--viral` implies `--all-genera`, facets on **family** rather than genus (the
output directory name is historical and wrong), disables the quality check, and
adds `--no-quality-check --no-check-files --batch-size 500` plus
`--max-missing-fraction 0.03` to each `p3x-create-blast-db`. That last one
exists because makeblastdb deterministically rejects a few records in some
viral sets; see "Verifying" below.

## Stage 3 — curated reference databases

Six databases in `ref/`: `bacteria-archaea` and `viral-reference`, each in the
three types. These are built by calling `p3x-create-blast-db` directly, with
`--reference --representative` to restrict to reference and representative
genomes rather than everything in the lineage.

The viral one is recorded verbatim in `ImplementationNotes.md`:

```bash
p3x-create-blast-db --no-quality-check --title "Viral Refs" --parallel 2 \
    --no-check-files --batch-size 500 \
    --reference --representative --taxon 10239 -f --viral \
    aa features $NEW/ref/viral-reference
```

> The `bacteria-archaea` invocation was not recorded anywhere in the tree and
> is **reconstructed** from the artifacts — taxa 2 (Bacteria) and 2157
> (Archaea), reference plus representative, which yields the 6,352 genomes and
> 6,232 taxa in the shipped `.glist`/`.taxlist`. Verify the genome count against
> the previous generation before trusting a rebuild:
>
> ```bash
> p3x-create-blast-db --title "Bacteria and Archaea" --parallel 2 \
>     --reference --representative --taxon 2 --taxon 2157 -f \
>     aa features $NEW/ref/bacteria-archaea
> ```

Repeat each for `dna features` and `dna contigs`.

**`--curated-directory` must match.** Stage 4 passes `--curated-directory ref`,
and the loader tests `$curated{dirname($dir)}` against the path *relative to the
database directory*. The argument therefore has to be exactly that relative
name. In `data.2022-0916` it did not match and all 4,058 `GenomeGroup` rows came
out `curated=0`, which makes the `NOT g.curated` filter in
`BlastDatabasesSQL::search_taxa` inert. Check it after stage 4 rather than
assuming.

## Stage 4 — the sqlite catalog

Everything here lives in this module's `scripts/`. The `p3x-` tools are wrapped
into `$KB_TOP/bin` and so are on PATH; `db-schema.sql` is data, read by path.

**The order is load-bearing:** `p3x-create-taxonomy-lineage` reads
`TaxonInDatabase`, which does not exist until `p3x-create-databases-lookup` has
walked the database tree.

**1. Download the NCBI taxdump and build the node table.**

    cd $NEW && p3x-create-taxonomy-nodes

Creates `taxdump-<YYYY-MM-DD>/` under the current directory containing
`nodes.tsv`, columns `tax_id,parent,rank`. The file is **comma**-separated
despite the name, which is why the import below needs `.mode csv`. The script
refuses to run if the dated directory already exists.

Do this at catalog-build time rather than reusing a taxdump fetched at the start
of the rebuild, so the taxonomy is never older than the genome set.

**2. Create the database and load the taxonomy.**

    sqlite3 $NEW/db.sqlite
    sqlite> .read scripts/db-schema.sql
    sqlite> .mode csv
    sqlite> .import taxdump-<YYYY-MM-DD>/nodes.tsv TaxNode

**3. Load the blast database catalog.**

    p3x-create-databases-lookup --curated-directory ref \
        --sqlite $NEW/db.sqlite $NEW blast.db

**4. Build the lineage lookup.**

    p3x-create-taxonomy-lineage $NEW/db.sqlite

`sqlite3` is not necessarily on PATH on the build host; `perl -MDBI` works for
inspection. For scale, `data.2022-0916` came out at: `BlastDatabase` 12,174 ·
`GenomeGroup` 4,058 · `GenomeInDatabase` 9,201,442 · `TaxNode` 2,443,146 ·
`TaxonInDatabase` 740,576 · `TaxonLineage` 2,919,811.

## Verifying before you publish

`TaxNode` should have exactly as many rows as `nodes.tsv` has lines.

`TaxonLineage` should have **no NULL `lineage_id`**. `p3x-create-taxonomy-lineage`
inserts a `(taxon_id, NULL)` row for any taxon in `TaxonInDatabase` that is
absent from `TaxNode` — the parent lookup returns undef, the `$p eq $cur`
termination guard is false, and one bogus row goes in before the walk ends.
Those taxa then match nothing and are silently unsearchable.
`data.2022-0916` has 63 of them, which is the symptom of a taxdump older than
the genome set.

For the databases themselves, three detectors. Note that **mean sequence length
does not work** — 1.5 Gaa of chimera is invisible diluted across 21.7M
sequences.

1. **Database records versus `.taxids` lines.** Assumption-free build-fidelity
   test, and as of the current tree `p3x-create-blast-db` performs it itself
   after makeblastdb: it reads `number-of-sequences` from the `.pjs`/`.njs`
   metadata, compares against the line count of `.taxids`, and on a mismatch
   removes `.taxids` so neither the skip guard nor the catalog loader will
   accept the database. The default tolerance is exact equality, because that is
   what clean bacterial builds produce; viral passes `0.03`. Reading the JSON
   sidecar is orders of magnitude cheaper than `blastdbcmd -info`.

2. **Features-database taxid coverage versus the same genus's contigs
   database.** Same genome list, so contigs clean plus features short means the
   intermediate fasta was truncated before makeblastdb ran. This is the detector
   that catches the case detector 1 cannot — `Acinetobacter.features.faa` in
   `data.2022-0916` is damaged yet its record count matches its map exactly.

3. **Longest protein over 50,000 aa, cross-checked against `fna == 3*faa (+3)`.**
   Consistency between the two means the giant record exists in both source
   files and is an upstream gene-calling artifact, not build corruption — 16
   genera look like this legitimately. Only the inconsistent ones are chimeras.

Two signals that look damning and are not: **`faa` records exceeding `fna`
records** occurs in the clean 2021 generation too and reflects the download
files, not the build; and **low taxid ratios on tiny genera** are small
denominators, where all three database types drop uniformly.

## Recovering a partial build

A generation that died partway is usually repairable in place. Restarting 1.8 TB
by reflex costs days; triage first, because the failure modes here are loud and
the damage is normally confined to the databases that actually failed.

Worked example, `blast-2026-0826`. The bacterial stage ran 08-26 to 08-29 against
a `P3DataAPI` that still paged with `start`/`rows`; deep paging without a total
ordering returns overlapping pages, so some genus enumerations came back with
**duplicate** genome ids. That defect is fixed (p3_core PR #24, cursor paging),
but the question was what the already-built 3.0 TB was worth.

**The bug was loud, not silent.** `makeblastdb -parse_seqids` rejects a duplicate
seqid, so an affected genus failed its build outright rather than quietly
shipping a wrong database. Triage bore that out:

    # 1. Which databases have sidecars but no blast index?
    cd $NEW/by-genus-bacterial
    for f in *.taxids; do b=${f%.taxids}
      [ -s "$b.psq" ] || [ -s "$b.nsq" ] || echo "$b"
    done

    # 2. Which genome lists contain duplicates?  (local, no network)
    for f in *.glist; do
      [ "$(wc -l < "$f")" = "$(sort -u "$f" | wc -l)" ] || echo "DUP $f"
    done

    # 3. Did any genus silently LOSE genomes?  (re-enumerate and diff)
    p3x-create-blast-db --dump --taxon <tax> > new.json
    # compare genome_ids against $NEW/by-genus-bacterial/<Genus>.features.faa.glist

Result: 14 of 11,125 databases lacked an index, 6 of 11,125 `.glist` files held
duplicates, and **every duplicated list was one of the failures** -- Campylobacter,
Escherichia, Salmonella, Staphylococcus, the largest genera, which is where deep
offsets hurt. Step 3 on the four largest genera found **zero lost genome ids**;
Salmonella's `+106` was a week of new submissions, not loss. So 11,111 databases
were sound and only the 14 needed rebuilding.

**Such a repair usually needs no `--overwrite`.** `makeblastdb` aborts before it
writes `.psq`/`.nsq`, and those are exactly the suffixes the skip guard tests
(`%blastdb_suffixes`), leaving only auxiliary `.pos`/`.pot`/`.ptf`/`.pto` behind.
So a plain rerun of the same three passes into the same directory skips the
11,111 good databases and rebuilds precisely the broken ones. Reach for
`--overwrite` only when the files are complete but the *content* is wrong -- the
detector-2 case, where a database and its taxid map agree with each other and
are both wrong.

Note what detector 1 cannot do here: `.taxids` is written from the same
enumeration as the database, so a count match proves the write was not truncated,
**not** that the genome list was complete. Only step 3, an independent
re-enumeration, tests that.

## Publishing

Copy or move the generation to `/vol/blastdb/bvbrc-service/data.<YYYY-MMDD>/`,
then update **both** stanzas in `deploy.cfg` — there are two,
`[homology_service]` and `[HomologyService]`, and each sets
`blast-db-search-path` and `blast-sqlite-db`. Missing one leaves half the
service on the old generation.

Keep the previous generation until the new one has served real queries.

## When the disk fills

Out of space does not stop the build cleanly. The forked children in
`process_from_download_files` write their fasta and taxid streams, the parent
`cat`s whatever is on disk into makeblastdb, and you get a short database with a
chimeric record spliced across each child-file boundary. Those writes are now
checked and the record-count gate in detector 1 catches the result, but the
build will still burn hours producing databases it then refuses — so check the
space first.

A related trap: `LPTScheduler::run` ignores the child exit status that
`Proc::ParallelLoop` records in `@Exit_Status`, so a child that dies is not by
itself visible to the parent. That is why the verification is a parent-side gate
on the finished artifact rather than only an error check in the child.
