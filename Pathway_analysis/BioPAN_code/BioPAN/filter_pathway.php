<?php
# Sanitization of form values to prevent injections and other attacks
$path = filter_input(INPUT_GET, 'path', FILTER_SANITIZE_STRING);
$selectedItems = filter_input(INPUT_GET, 'selectedItems', FILTER_SANITIZE_STRING);
$userDir = filter_input(INPUT_GET, 'userDir', FILTER_SANITIZE_STRING);
$outFilePath = "$userDir/biopan/";

# Write selected lipids to csv file
$file = "$userDir/biopan/selected_items.csv";
$fp = fopen($file, 'w');
fwrite($fp, $selectedItems . PHP_EOL);
fclose($fp);

# Execute R script
exec("/lipidmaps/R-4.0.2/bin/Rscript R/filter_pathway.r $outFilePath $path &>> $userDir/error.log");
?>
