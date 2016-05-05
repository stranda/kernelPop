TEX2DVI	= texi2dvi --pdf

all: kernelPop-intro.pdf kernelPop-intro.R

kernelPop-intro.pdf: kernelPop-intro.Rnw
	echo 'Sweave("kernelPop-intro.Rnw")'|R --no-save
	${TEX2DVI} kernelPop-intro.tex
	cp kernelPop-intro.pdf ../
	rm kernelPop-intro-*.eps kernelPop-intro-*.pdf
kernelPop-intro.R: kernelPop-intro.Rnw
	echo 'Stangle("kernelPop-intro.Rnw")'|R --no-save
	cp kernelPop-intro.R ../
