
ZIP_FILES=$(wildcard ./zips/*)
ZIP300=$(wildcard ./zips/*300*)
ZIP800=$(wildcard ./zips/*800*)
# ZIP_BASE=$(subst ./zips/,,$(ZIP_FILES))
ZIP_BASE=$(patsubst ./zips/%-300-500.zip,%,$(ZIP300))

DATA300=$(patsubst %,exp1/%-300.csv,$(ZIP_BASE))
DATA800=$(patsubst %,exp1/%-800.csv,$(ZIP_BASE))
ADD300=$(patsubst %,exp2/%-300.csv,$(ZIP_BASE))
ADD800=$(patsubst %,exp2/%-800.csv,$(ZIP_BASE))

EXP3=$(patsubst exp2/%,exp3/%, $(ADD300) $(ADD800))
EXP4=$(patsubst exp2/%,exp4/%, $(ADD300) $(ADD800))
$(info ZIP_FILES Is [${ZIP_FILES}])

$(info DATA300 Is [${DATA300}])
$(info ZIP_BASE Is [${ZIP_BASE}])
$(info EXP3 Is [${EXP3}])


.PHONY: 

.ONESHELL:
unzips: ${ZIP_FILES}
	mkdir unzips
	mkdir unzips/n300
	mkdir unzips/n800
	for F in ${ZIP300}; do 
		unzip -j $${F}  -d unzips/n300 
	done
	for F in ${ZIP800}; do 
		unzip -j $${F}  -d unzips/n800 
	done

data:
		mkdir data

data/separation.csv: data separation.R
		Rscript -e "source('separation.R'); df <- map(files, check_sep) |> list_rbind(); write.csv(df, '$@', quote=F, row.names = F)"

data/multi-acc.csv: data lda.R
		Rscript -e "source('lda.R'); write_multi(1, tol = 1/2^(4:12))"

data/triples.csv: data lda.R
		Rscript -e "source('lda.R'); write_triples(1)"

tab-sep.tex: data/separation.csv
		Rscript -e "source('tab-sep.R')"

tab-multi.tex: data/multi-acc.csv
		Rscript -e "source('tab-multi.R')"

data/ex2-%-300.csv: unzips data

exp1/%-300.csv: exp1
ifeq (,$(wildcard $@))
				Rscript -e "source('msvm.R');df<-sah1(300,'$*');write.csv(df, '$@', quote=F, row.names=F)"
else
				echo "Missing rule"
endif

paper1.pdf: paper1.tex tab-sep.tex tab-multi.tex
		pdflatex paper1.tex

graph1.pdf graph2.pdf: graph12.R
		Rscript -e "source('graph12.R')"

all: $(DATA300) $(DATA800) $(ADD300) $(ADD800) $(EXP3) $(EXP4)

clean:
		rmdir -rf unzips
		rm tab*.tex
		rmdir -rf exp1
		rmdir -rf exp2
