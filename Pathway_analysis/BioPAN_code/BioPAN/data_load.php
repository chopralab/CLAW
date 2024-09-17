<?php
$userID = filter_input(INPUT_POST, "id", FILTER_SANITIZE_STRING);
$userID = preg_replace("[^A-Za-z0-9]", "", $userID);
# User's root folder
$userDir = "/lipidmaps/temp/$userID";
if (empty($userID)) {
    include("includes/top.php");
    
    echo "<center><h1>BioPAN: Bioinformatics Methodology For Pathway Analysis</h1></center><br><br><br>
          <strong>ERROR:</strong> Unrecognised user \"$userID\".<br><br>
          Please go back to the <a href='index.php'>index</a> to get a new identifier.";
    include("includes/bottom.php");
    exit;
}
# Change the permissions for newly created files
umask(0022);
# If the user's root folder already exists, delete it to start a new session
if (file_exists($userDir)) {
    # Recursive function to delete every subfolder and file in a given directory
    function rrmdir($dir) {
        if (is_dir($dir)) {
            $objects = scandir($dir);
            foreach ($objects as $obj) {
                if (($obj != '.') && ($obj != '..')) {
                    if (filetype("$dir/$obj") == 'dir') {
                        rrmdir("$dir/$obj");
                    } else {
                        unlink("$dir/$obj");
                    }
                }
            }
            reset($objects);
            rmdir($dir);
        }
    }
    rrmdir($userDir);
}
# Set a time limit to avoid getting stuck when uploading really large files
set_time_limit();
# When loading a new file, every user's previous folder has been deleted
mkdir($userDir);
mkdir("$userDir/input");
mkdir("$userDir/config");

# Set the file path where the data CSV file will be saved
$dstFilePath = "$userDir/input/input.csv";

# Sanitization of form values to prevent injections and other attacks
$dataType = filter_input(INPUT_POST, "dataType", FILTER_SANITIZE_STRING);

$file = False;

if($dataType == "demo") {
    $demoFilePath = "resources/sample.csv";
    copy($demoFilePath, $dstFilePath);
    $params = array("srcFileName" => "demo sample", "lipidlynxx" => "no");
    $file = True;
} elseif($dataType == "smallDemo") {
    $demoFilePath = "resources/smallSample.csv";
    copy($demoFilePath, $dstFilePath);
    $params = array("srcFileName" => "Short demo sample", "lipidlynxx" => "no");
    $file = True;

} elseif($dataType == 'local') {
    $tmpFilePath = $_FILES['dataFile']['tmp_name'];
    $srcFileName = $_FILES['dataFile']['name'];
    if (isset($tmpFilePath) && move_uploaded_file($tmpFilePath, $dstFilePath)) {
        $params = array("srcFileName" => $srcFileName, "lipidlynxx" => "yes");
        $file = True;
    } else {
        # Something has gone wrong whilst uploading or moving the data file
        include("includes_LM/top.php");
        echo "<center><h1>BioPAN: Bioinformatics Methodology For Pathway Analysis</h1></center><br><br><br>
              <div><strong>There was an error uploading the data file <font color=#ff000>
              $srcFileName</font></strong>.<br><br>Please go back and try again!</div>";
        include("includes_LM/bottom.php");
    }
} else {
    echo "Unexpected value in form. Please, go back and try again.";
}

# The file has been successfully uploaded and moved into the user's input folder
if($file = TRUE){
    # Save the source file name in a JSON file so they will be still available if the user leaves and comes back
    $fp = fopen("$userDir/config/biopan_params.json", 'w');
    fwrite($fp, json_encode($params));
    fclose($fp);
    # Go to the next step of the process
    header("Location: parse_data.php?id=$userID");
}

exit;
?>
