<?php
include("$root/includes/top.php");
$userID=filter_input(INPUT_GET, 'id', FILTER_SANITIZE_STRING);
$userID = preg_replace("[^A-Za-z0-9]", "", $userID);
echo "<center><h1>LipidFinder: LC/MS Analysis Workflow</h1></center><br><br><br>
          <strong>ERROR:</strong> Unrecognised user \"$userID\".<br><br>
          Please go back to the <a href='index.php'>index</a> to get a new identifier.";
include("$root/includes/bottom.php");
?>
