<?php
# Sanitization of form values to prevent injections and other attacks
$path = filter_input(INPUT_GET, 'path', FILTER_DEFAULT, FILTER_REQUIRE_ARRAY);

$pathArray = array();

for($i = 0; $i < count($path); $i++){
    $string = file_get_contents($path[$i]);
    $json_obj = json_decode($string, true);
    array_push($pathArray, json_encode($json_obj));
}

header('Content-Type: application/json');
echo json_encode($pathArray);
?>
