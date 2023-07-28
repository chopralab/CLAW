<?php
# Sanitization of form values to prevent injections and other attacks
$path = filter_input(INPUT_GET, 'path', FILTER_SANITIZE_STRING);
$string = file_get_contents($path);
echo $string;
?>
