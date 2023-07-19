<?php
$userID = filter_input(INPUT_GET, "id", FILTER_SANITIZE_STRING);
$userID = preg_replace("[^A-Za-z0-9]", "", $userID);
# User's root folder
$userDir = "/lipidmaps/temp/$userID";
# Show error if the page has been accessed with a wrong ID
if (empty($userID) || !file_exists($userDir)) {
    include("includes_LM/top.php");
    echo "<center><h1>BioPAN: BIOPAN's Workflow </h1></center><br><br><br>
          <strong>ERROR:</strong> Unrecognised user \"$userID\".<br><br>
          Please go back to the <a href='index.php'>index</a> to get a new identifier.";
    include("includes_LM/bottom.php");
    exit;
}
# Change the permissions for newly created files
umask(0022);
if (!file_exists("$userDir/biopan")) {
    # Create the output folder if the user's "workspace" has been reset
    mkdir("$userDir/biopan");
} else {
    # Remove the existing files to keep track of the new process
    $objects = scandir("$userDir/biopan");
    foreach ($objects as $obj) {
        if (($obj != '.') && ($obj != '..') && is_file("$userDir/biopan/$obj")) {
            unlink("$userDir/biopan/$obj");
        }
    }
}

mkdir("$userDir/biopan/reaction");
$inputPath = "$userDir/input/input.csv";
$inputPathClean = "$userDir/input/input_clean.csv";
$outputPath = "$userDir/biopan/";
$outputLipidLynxX = "$userDir/input/input_LipidLynxX.csv";

if (file_exists($outPathwayFilePath)) {
    unlink($outPathwayFilePath);
}
if (file_exists($outNetFilePath)) {
    unlink($outNetFilePath);
}

# Retrieve information on running LipidLynxX or not
$str = file_get_contents("$userDir/config/biopan_params.json");
$privateParams = json_decode($str, TRUE);

# Load information on running LipidLynxX or not (don't need to run it for the demo files)
$str = file_get_contents("$userDir/config/biopan_params.json");
$privateParams = json_decode($str, TRUE);
$lipidlynxx = $privateParams['lipidlynxx'];

# Execute LipidLynxX
if($lipidlynxx == "yes"){
    shell_exec("./runLipidLynxX.sh $inputPath $outputLipidLynxX &>> $userDir/error1.log");
}

# Execute R script
exec("/lipidmaps/R-4.0.2/bin/Rscript R/parse_data.r $inputPath $inputPathClean $outputPath $lipidlynxx $outputLipidLynxX &>> $userDir/error.log");
header("Location: summary.php?id=$userID");
exit;
?>
