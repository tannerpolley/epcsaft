# Build: cd docs/scientific/latex && latexmk   (XeLaTeX; main.pdf and supplement.pdf in builds/, cross-referenced through xr-hyper; pdfLaTeX embeds bitmap sans fonts)
@default_files = ("main.tex", "supplement.tex");
$pdf_mode = 5;
$out_dir = 'builds';
$aux_dir = 'builds';
$ENV{'TEXINPUTS'} = './styles//;' . ($ENV{'TEXINPUTS'} || '');
