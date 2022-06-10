TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime
SERVER_SPEC = HomologyService.spec

APP_SERVICE = app_service

#
# This layout is specfic to the current (January 2016) layout of 
# the NCBI BLAST distribution.
#
BLAST_VERSION = 2.13.0
BLAST_BASE = ncbi-blast-$(BLAST_VERSION)+
BLAST_FTP_SRC = ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$(BLAST_VERSION)/$(BLAST_BASE)-x64-linux.tar.gz
BLAST_FTP_FILE = $(notdir $(BLAST_FTP_SRC))
BLAST_DEPLOY_DIR = $(TARGET)/services/$(SERVICE)/bin
BLAST_DB_SEARCH_PATH = /vol/blastdb/bvbrc-service

SERVICE_MODULE = lib/Bio/KBase/HomologyService/Service.pm

SERVICE = homology_service
SERVICE_NAME = HomologyService
SERVICE_PORT = 7134

#ASYNC_SERVICE_PORT = 7135
ASYNC_SERVICE_PSGI = HomologyServiceAsync.psgi

SERVICE_URL = https://kbase.us/services/$(SERVICE)

SERVICE_NAME = HomologyService
SERVICE_NAME_PY = $(SERVICE_NAME)

SERVICE_PSGI_FILE = $(SERVICE_NAME).psgi

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))


ifdef TEMPDIR
TPAGE_TEMPDIR = --define kb_tempdir=$(TEMPDIR)
endif

TPAGE_ARGS = --define kb_top=$(TARGET) \
	--define kb_runtime=$(DEPLOY_RUNTIME) \
	--define kb_service=$(SERVICE) \
	--define kb_service_name=$(SERVICE_NAME) \
	--define kb_service_port=$(SERVICE_PORT) \
	--define kb_psgi=$(SERVICE_PSGI_FILE) \
	--define kb_async_service_port=$(ASYNC_SERVICE_PORT) \
	--define kb_async_psgi=$(ASYNC_SERVICE_PSGI) \
	--define blast_db_search_path=$(BLAST_DB_SEARCH_PATH) \
	--define blast_sqlite_db=$(BLAST_SQLITE_DB) \
	$(TPAGE_TEMPDIR)

TESTS = $(wildcard t/client-tests/*.t)

all: build-libs build-blast bin compile-typespec service

test:
	# run each test
	echo "RUNTIME=$(DEPLOY_RUNTIME)\n"
	for t in $(TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

service:

compile-typespec: Makefile
	mkdir -p lib/biop3/$(SERVICE_NAME_PY)
	touch lib/biop3/__init__.py #do not include code in biop3/__init__.py
	touch lib/biop3/$(SERVICE_NAME_PY)/__init__.py 
	mkdir -p lib/javascript/$(SERVICE_NAME)
	compile_typespec \
--patric \
		--psgi $(SERVICE_PSGI_FILE) \
		--impl Bio::KBase::$(SERVICE_NAME)::%sImpl \
		--service Bio::KBase::$(SERVICE_NAME)::Service \
		--client Bio::KBase::$(SERVICE_NAME)::Client \
		--py biop3/$(SERVICE_NAME_PY)/client \
		--js javascript/$(SERVICE_NAME)/Client \
		--url $(SERVICE_URL) \
		--enable-retries \
		$(SERVER_SPEC) lib
	-rm -f lib/$(SERVER_MODULE)Server.py
	-rm -f lib/$(SERVER_MODULE)Impl.py
	-rm -f lib/CDMI_EntityAPIImpl.py

bin: $(BIN_PERL) $(BIN_SERVICE_PERL)

deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service
deploy-client: build-libs compile-typespec deploy-docs deploy-libs deploy-scripts 

build-libs:
	$(TPAGE) $(TPAGE_BUILD_ARGS) $(TPAGE_ARGS) Config.pm.tt > lib/Bio/P3/HomologySearch/Config.pm

deploy-service: deploy-dir deploy-libs deploy-service-scripts deploy-blast deploy-specs
	$(TPAGE) $(TPAGE_ARGS) service/start_service.tt > $(TARGET)/services/$(SERVICE)/start_service
	chmod +x $(TARGET)/services/$(SERVICE)/start_service
	$(TPAGE) $(TPAGE_ARGS) service/stop_service.tt > $(TARGET)/services/$(SERVICE)/stop_service
	chmod +x $(TARGET)/services/$(SERVICE)/stop_service

deploy-specs:
	mkdir -p $(TARGET)/services/$(APP_SERVICE)
	rsync -arv app_specs $(TARGET)/services/$(APP_SERVICE)/.

deploy-service-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib ; \
	export PATH_PREFIX=$(TARGET)/services/$(SERVICE)/bin:$(TARGET)/services/cdmi_api/bin; \
	for src in $(SRC_SERVICE_PERL) ; do \
	        basefile=`basename $$src`; \
	        base=`basename $$src .pl`; \
	        echo install $$src $$base ; \
	        cp $$src $(TARGET)/plbin ; \
	        $(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/bin/$$base ; \
	done

#
# We use a captive blast build in order to tightly control versioning.
#

build-blast: blast.deploy/$(BLAST_BASE)/bin/blastp

deploy-blast: build-blast
	mkdir -p $(BLAST_DEPLOY_DIR)
	cp blast.deploy/$(BLAST_BASE)/bin/* $(BLAST_DEPLOY_DIR)

blast.deploy/$(BLAST_BASE)/bin/blastp:
	if [ ! -s $(BLAST_FTP_FILE) ] ; then \
		curl --fail -o $(BLAST_FTP_FILE) $(BLAST_FTP_SRC);  \
	fi
	rm -rf blast.deploy
	mkdir blast.deploy
	tar -C blast.deploy -x -v -f $(BLAST_FTP_FILE)
	rm -f blast.bin
	ln -s blast.deploy/$(BLAST_BASE)/bin blast.bin

doit:
	-ldso=$$KB_RUNTIME/glibc-2.17/lib/ld-linux-x86-64.so.2 ; \
	if grep -i 'release 6' /etc/centos-release && [[ -f $$ldso ]] ; then \
	    for f in blast.deploy/$(BLAST_BASE)/bin/* ; do \
		echo "Patch $$f"; \
		patchelf --set-interpreter $$ldso $$f; \
	    done \
	fi

deploy-monit:
	$(TPAGE) $(TPAGE_ARGS) service/process.$(SERVICE).tt > $(TARGET)/services/$(SERVICE)/process.$(SERVICE)

deploy-docs:
	-mkdir doc
	-mkdir $(SERVICE_DIR)
	-mkdir $(SERVICE_DIR)/webroot
	mkdir -p doc
	$(DEPLOY_RUNTIME)/bin/pod2html -t "Homology Search Service API" lib/Bio/KBase/HomologyService/HomologyServiceImpl.pm > doc/homologyservice_impl.html
	cp doc/*html $(SERVICE_DIR)/webroot/.

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/webroot ] ; then mkdir $(SERVICE_DIR)/webroot ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

$(BIN_DIR)/%: service-scripts/%.pl $(TOP_DIR)/user-env.sh
	$(WRAP_PERL_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/%: service-scripts/%.py $(TOP_DIR)/user-env.sh
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

include $(TOP_DIR)/tools/Makefile.common.rules
