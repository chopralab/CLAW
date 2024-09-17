<?php

function lmGetFileDownloadUrl($name, $ext) {
    return "/files/?file={$name}&ext={$ext}";
}