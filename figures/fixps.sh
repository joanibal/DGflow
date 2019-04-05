
NAME="machp2"

ps2pdf "$NAME.ps"
pdftk "$NAME.pdf" cat 1-endeast output "$NAME+out.pdf"
pdfcrop "$NAME+out.pdf" "$NAME.pdf"