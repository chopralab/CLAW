<?php
//Fetching values from url and storing it in PHP variables.
// Sanitization of form values to prevent injections and other attacks
$userID=$_POST['userID'];
$userID = preg_replace("[^A-Za-z0-9]", "", $userID);
$str_samples = $_POST['samples'];
$str_groups = $_POST['groups'];


# User's root folder
$userDir = "/lipidmaps/temp/$userID";
$response = "";
# Show error if the page has been accessed with a wrong ID
if (empty($userID) || !file_exists($userDir)) {
    $response = "user_error.php?id=$userID";
}else{
    $p_value = 0.05;
    $is_paired = 'FALSE';
    $pathway_level = "default";
    $inFilePath = "$userDir/input/input_clean.csv";
    $outFilePath = "$userDir/biopan/";
    $summaryFilePath = "$userDir/biopan/summary.json";
    # Execute R script
    exec("/lipidmaps/R-4.0.2/bin/Rscript R/process_data.r $inFilePath $outFilePath $summaryFilePath $str_samples $str_groups $p_value $is_paired $pathway_level&>> $userDir/error.log");
    $response = "pathway_analysis.php?id=$userID";
}
echo $response;
?>
