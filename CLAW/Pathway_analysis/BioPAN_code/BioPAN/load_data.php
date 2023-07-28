<?php
$userID = substr(str_shuffle("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"), 0, 16);
# User's root folder
$userDir = "/lipidmaps/temp/$userID";
include("includes_LM/top.php");
?>

<html>
<head>
    <title>LIPID MAPS® Lipidomics Gateway</title>
    <meta name="viewport" content="width=device-width, initial-scale=1  minimum-scale=1.0">
    <meta charset="utf-8">
    <meta http-equiv="x-ua-compatible" content="ie=edge">

    <script src='js/librairies/jquery.min.js'></script>
    <script src="js/upload.js"></script>
    <link rel="stylesheet" href="css/processing.css" />
    <link rel='stylesheet' href='https://fonts.googleapis.com/icon?family=Material+Icons'>
</head>

<body>
    <center><h1>BioPAN's Workflow</h1></center><br>

    <div id='lf_step1' style='line-height: 1.3;'>
        <h2>STEP 1: Load a data file</h2>
        
        <br>

        <!-- Loading figure -->
        <div class="overlay" style="display:none;">
            <div class="loading">
                <p>Uploading file...</p>
            </div>
        </div>

        <p>There are guidelines on how prepare your data <a href='doc/step1.html#how-to-prepare-your-data' target='_blank'>here</a>. <br>
        After loading a file, LipidLynxX<sup><b>1</b></sup> will be firslty launched in the background to convert the lipid molecular species naming used in the file to the BioPAN nomenclature.
        </p>    
        <br>

        <ul>
            <li style='padding-bottom:1rem;'><h5>Option 1: select a data file for analysis</h5></li>
            <form action='data_load.php' method='POST' enctype='multipart/form-data'>
                <input type='hidden' name='id' value='<?php echo "$userID"; ?>'>
                <input type='hidden' name='dataType' value='local'>

                <div id='tbl_user_data'>
                    <input type='file' id='dataFile' name='dataFile' onchange='return fileValidation()' required>
                    <input id="run" type='submit' value='Load data file' style='margin-top:15px; display: flex; align-items: center; justify-content: center; height: 26px; width: 140px; font-weight: bold; cursor: pointer;'>
                </div>
            </form>

        <!-- Launch BioPAN if a dataset is selected -->
        <script>
        document.getElementById('run').onclick = function() {
            if(document.getElementById('dataFile').value != ''){
                $(".overlay").show(); $(".loading").css("top","54%"); $(".loading").css("left","47%");
            }
        }
        </script>
        
        <br>
        <br>
        <br>

        <li style='padding-bottom:1rem;'><h5>Option 2: try BioPAN with a demonstration file</h5></li>
        <form action='data_load.php' method='POST' enctype='multipart/form-data'>
            <input type='hidden' name='id' value='<?php echo "$userID"; ?>'>
            <div id='tbl_demo_data'>
                <select name="dataType">
                    <option value="smallDemo">Small dataset</option>
                    <option value="demo">Complete dataset</option>
                </select>

                <input type='submit' value='Load demo file' onclick='$(".overlay").show(); $(".loading").css("top","54%"); $(".loading").css("left","47%");' style='margin-top:15px; display: flex; align-items: center; justify-content: center; height: 26px; width: 140px; font-weight: bold; cursor: pointer;'>
            </div>
        </form>

    </ul>

        <br>
        
        <p> The small demonstration file is downloadable <button onclick="download('resources/smallSample.csv','BioPAN_small_demo.csv');" style="text-decoration: underline; color: #1F5CAA; padding:0; border:none; background:none; cursor: pointer;">here</button> and the complete dataset<sup><b>2</b></sup> <button onclick="download('resources/sample.csv','BioPAN_demo.csv');" style="text-decoration: underline; color: #1F5CAA; padding:0; border:none; background:none; cursor: pointer;">here</button>.<p>
        
        <br>
        <br>
        <p><sup><b>1</b></sup><strong> LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets.</strong> Zhixu Ni, Maria Fedorova. <br>
            <a href='https://www.biorxiv.org/content/10.1101/2020.04.09.033894v2' target='_blank'>bioRxiv, 2020.04.09.033894.</a></p>
        <br>
        <p><sup><b>2</b></sup><strong> A nutritional memory effect counteracts the benefits of dietary restriction in old mice.</strong> Oliver Hahn, Lisa F. Drews, An Nguyen et al. <br>
            <a href='https://www.nature.com/articles/s42255-019-0121-0' target='_blank'>Nat Metab 1, 1059-1073 (2019).</a></p>
        <br>
    </div>
</body>
</html>

<?php
include("includes_LM/bottom.php");
?>
