# Build: cd docs/scientific/latex && latexmk main.tex   (XeLaTeX, PDF in builds/; pdfLaTeX embeds bitmap sans fonts)
$pdf_mode = 5;
$out_dir = 'builds';
$aux_dir = 'builds';
$ENV{'TEXINPUTS'} = './styles//;' . ($ENV{'TEXINPUTS'} || '');
