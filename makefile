
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
		unzip -n -j $${F}  -d unzips/n300 
	done
	for F in ${ZIP800}; do 
		unzip -n -j $${F}  -d unzips/n800 
	done

data: unzips
		mkdir -p data

graphs: 
		mkdir -p graphs


data/separation.csv: data separation.R
		Rscript -e "source('separation.R'); plan(multicore, workers = 4); df <- future_map(files, check_sep) |> list_rbind(); write.csv(df, '$@', quote=F, row.names = F)"

data/multi-acc.csv: data lda.R
		Rscript -e "source('lda.R'); write_multi(1, tol = 1/2^(4:12))"
		
data/exp4.csv: data lda.R  exp4.R
		Rscript -e "source('exp4.R')"

data/triples.csv: data lda.R
		Rscript -e "source('lda.R'); write_triples(1)"

data/exp2.csv: data exp2-summary.R
		Rscript -e "source('exp2-summary.R')"


data/bc3_%.csv: data lda.R
		Rscript -e "source('lda.R');pt<-par_triples(limit = 1000, score =
		'$*');print(class(pt));print(class(pt\$$df));write.csv(pt\$$df, '$@', quote=F, row.names=F)"
 
exp2-detail.pdf: exp2-detail.R graphs theme.R
		Rscript -e "source('exp2-detail.R')"

exp2-plot2.pdf: exp2-plot2.R graphs theme.R data/exp2.csv
		Rscript -e "source('exp2-plot2.R')"

exp1-plot1.pdf: exp1-plot1.R data/multi-acc.csv theme.R
		Rscript -e "source('exp1-plot1.R')"

exp4-truth.pdf: data/exp4.csv exp4-truth.R
		Rscript -e "source('exp4-truth.R')"

tab-sep.tex: data/separation.csv tab-sep.R
		Rscript -e "source('tab-sep.R')"

tab-step2.tex: data/triples.csv tab-step2.R
		Rscript -e "source('tab-step2.R')"

tab-multi.tex: data/multi-acc.csv
		Rscript -e "source('tab-multi.R')"

tab-exp4.tex: data/exp4.csv tab-exp4.R
		Rscript -e "source('tab-exp4.R')"

data/ex2-%-300.csv: unzips data

exp1/%-300.csv: exp1
ifeq (,$(wildcard $@))
				Rscript -e "source('msvm.R');df<-sah1(300,'$*');write.csv(df, '$@', quote=F, row.names=F)"
else
				echo "Missing rule"
endif

TEX_FILES= intro.tex theory.tex experiment.tex coupling.tex discussion.tex appendix.tex tab-sep.tex tab-step2.tex tab-exp4.tex
FIG_FILES= exp2-plot2.pdf  exp4-truth.pdf

paper.pdf: paper.tex tab-sep.tex tab-step2.tex exp2-plot2.pdf exp2-summary.pdf exp2-detail.pdf exp1-plot1.pdf exp4-truth.pdf paper.bib tab-exp4.tex 
		pdflatex paper.tex
		bibtex paper
		pdflatex paper.tex
		pdflatex paper.tex


sp1.pdf: springer.tex $(TEX_FILES) $(FIG_FILES)
		pdflatex "\\def\\mysecret{1} \\input{springer.tex}" 
		mv springer.pdf sp1.pdf

sp2.pdf: springer.tex $(TEX_FILES) paper.bib $(FIG_FILES)
		pdflatex "\\def\\mysecret{2} \\input{springer.tex}" 
		bibtex springer
		pdflatex "\\def\\mysecret{2} \\input{springer.tex}" 
		pdflatex "\\def\\mysecret{2} \\input{springer.tex}" 
		mv springer.pdf sp2.pdf

sp3.pdf: springer.tex  $(TEX_FILES) paper.bib $(FIG_FILES)
		pdflatex "\\def\\mysecret{3} \\input{springer.tex}" 
		bibtex springer
		pdflatex "\\def\\mysecret{3} \\input{springer.tex}" 
		pdflatex "\\def\\mysecret{3} \\input{springer.tex}" 
		mv springer.pdf sp3.pdf


all: $(DATA300) $(DATA800) $(ADD300) $(ADD800) $(EXP3) $(EXP4)

clean:
		rm -rf unzips
		rm -f tab*.tex *.toc *.out *.pdf *.aux *.log
