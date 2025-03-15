
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


exp1: 
		mkdir exp1

exp2: 
		mkdir exp2

exp%:
		mkdir $@

data:
		mkdir data

data/separation.csv: data separation.R
		Rscript -e "source('separation.R'); df <- map(files, check_sep) |> list_rbind(); write.csv(df, '$@', quote=F, row.names = F)"

data/multi.csv: data lda.R
		Rscript -e "source('lda.R'); write_multi(1)"

data/ex2-%-300.csv: unzips data

exp1/%-300.csv: exp1
ifeq (,$(wildcard $@))
				Rscript -e "source('msvm.R');df<-sah1(300,'$*');write.csv(df, '$@', quote=F, row.names=F)"
else
				echo "Missing rule"
endif


exp1/%-800.csv: exp1
ifneq (,$(wildcard $@))
		Rscript -e "source('msvm.R');df<-sah1(800,'$*');write.csv(df, '$@', quote=F, row.names=F)"
endif

exp2/%-300.csv: exp2
ifneq (,$(wildcard $@))
		Rscript -e "source('msvm.R');df<-sah1(300,'$*',lapl=-2:4);write.csv(df, '$@', quote=F, row.names=F)"
endif

exp3/%-300.csv: exp3
		Rscript -e "source('msvm.R');df<-sah1(300,'$*',lapl=-3:6);write.csv(df, '$@', quote=F, row.names=F)"


exp2/%-800.csv: exp2
ifneq (,$(wildcard $@))
		Rscript -e "source('msvm.R');df<-sah1(800,'$*', lapl=-2:4);write.csv(df, '$@', quote=F, row.names=F)"
endif

exp3/%-800.csv: exp3
		Rscript -e "source('msvm.R');df<-sah1(800,'$*', lapl=-3:6);write.csv(df, '$@', quote=F, row.names=F)"

exp4/%-800.csv: exp4
		Rscript -e "source('msvm.R');df<-sah1(800,'$*', cent = c(T,F));write.csv(df, '$@', quote=F, row.names=F)"

exp4/%-300.csv: exp4
		Rscript -e "source('msvm.R');df<-sah1(300,'$*', cent = c(T,F));write.csv(df, '$@', quote=F, row.names=F)"


paper1.pdf: paper1.tex exp1anova.tex graph1.pdf graph2.pdf
		pdflatex paper1.tex

graph1.pdf graph2.pdf: graph12.R
		Rscript -e "source('graph12.R')"

exp1anova.tex: exp1.R
		Rscript -e "source('exp1.R')"

exp2anova.tex: exp1.R
		Rscript -e "source('exp2.R')"


all: $(DATA300) $(DATA800) $(ADD300) $(ADD800) $(EXP3) $(EXP4)

clean:
		rmdir -rf unzips
		rmdir -rf exp1
		rmdir -rf exp2
