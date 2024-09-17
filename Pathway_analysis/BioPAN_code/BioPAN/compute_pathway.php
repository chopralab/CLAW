<?php
# Sanitization of form values to prevent injections and other attacks
$pwLevel= filter_input(INPUT_GET, 'pw_level', FILTER_SANITIZE_STRING);
$pValue= filter_input(INPUT_GET, 'p_value', FILTER_SANITIZE_STRING);
$isPaired= filter_input(INPUT_GET, 'is_paired', FILTER_SANITIZE_STRING);
$userDir = filter_input(INPUT_GET, 'userDir', FILTER_SANITIZE_STRING);
$outFilePath = "$userDir/biopan/";

# Execute R script
exec("/lipidmaps/R-4.0.2/bin/Rscript R/compute_pathway.r $outFilePath $pwLevel $pValue $isPaired &>> $userDir/error.log");
?>
