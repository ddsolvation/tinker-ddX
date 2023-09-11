TARGETS = dmso_cosmo dmso_pcm dmso_pcm_einf dmso_lpb dmso_lpb_einf dmso_lpb_kinf dmso_lpb_kzero
LIMITS = dmso_pcm_einf_limit dmso_lpb_einf_limit dmso_lpb_kinf_limit dmso_lpb_kzero_limit

XYZ_FILES := $(addsuffix .xyz, $(TARGETS))
KEY_FILES := $(addsuffix .key, $(TARGETS))
LOG_FILES := $(addsuffix .log, $(TARGETS))

.PHONY: clean all $(TARGETS) $(LIMITS)

all: clean $(TARGETS) $(LIMITS)

dmso_pcm_einf_limit: dmso_cosmo.log dmso_pcm_einf.log
	./compare.py dmso_cosmo.log dmso_pcm_einf.log

dmso_lpb_einf_limit: dmso_lpb_einf.log dmso_cosmo.log
	./compare.py dmso_lpb_einf.log dmso_cosmo.log

dmso_lpb_kinf_limit: dmso_lpb_kinf.log dmso_cosmo.log
	./compare.py dmso_lpb_kinf.log dmso_cosmo.log loose

dmso_lpb_kzero_limit: dmso_lpb_kzero.log dmso_pcm.log
	./compare.py dmso_lpb_kzero.log dmso_pcm.log loose

$(TARGETS): % : %.log %.ref
	./compare.py $@.log $@.ref

$(LOG_FILES): %.log: %.key %.xyz
	./analyze.x $(basename $@) E > $@

$(XYZ_FILES): %.xyz: dmso.xyz
	cp $< $@

dmso_cosmo.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddcosmo/g' $@
	sed -i 's/EPS/1e10/g' $@
	sed -i 's/KAPPA/0.104/g' $@
	sed -i 's/FMM/0/g' $@

dmso_pcm.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddpcm/g' $@
	sed -i 's/EPS/80.0/g' $@
	sed -i 's/KAPPA/0.104/g' $@
	sed -i 's/FMM/0/g' $@

dmso_pcm_einf.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddpcm/g' $@
	sed -i 's/EPS/1e10/g' $@
	sed -i 's/KAPPA/0.104/g' $@
	sed -i 's/FMM/0/g' $@

dmso_lpb.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddlpb/g' $@
	sed -i 's/EPS/80.0/g' $@
	sed -i 's/KAPPA/0.104/g' $@
	sed -i 's/FMM/0/g' $@

dmso_lpb_einf.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddlpb/g' $@
	sed -i 's/EPS/1e10/g' $@
	sed -i 's/KAPPA/0.104/g' $@
	sed -i 's/FMM/0/g' $@

dmso_lpb_kinf.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddlpb/g' $@
	sed -i 's/EPS/80.0/g' $@
	sed -i 's/KAPPA/10.0/g' $@
	sed -i 's/FMM/0/g' $@

dmso_lpb_kzero.key: dmso.key
	cp dmso.key $@
	sed -i 's/METHOD/ddlpb/g' $@
	sed -i 's/EPS/80.0/g' $@
	sed -i 's/KAPPA/1e-6/g' $@
	sed -i 's/FMM/0/g' $@

clean:
	rm -f *.log $(XYZ_FILES) $(KEY_FILES)
