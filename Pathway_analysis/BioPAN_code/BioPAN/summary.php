<?php
include("includes_LM/apache_request_headers.php");

$userID = filter_input(INPUT_GET, 'id', FILTER_SANITIZE_STRING);
$userID = preg_replace("[^A-Za-z0-9]", "", $userID);
# User's root folder
$userDir = "/lipidmaps/temp/$userID";
# Show error if the page has been accessed with a wrong ID
if (empty($userID) || !file_exists($userDir)) {
    include("includes_LM/top.php");
    echo "<center><h1>BioPAN: BioPAN's Workflow </h1></center><br><br><br>
          <strong>ERROR:</strong> Unrecognised user \"$userID\".<br><br>
          Please go back to the <a href='index.php'>index</a> to get a new identifier.";
    include("includes_LM/bottom.php");
    exit;
}

# Load data file name saved in data_load.php
$str = file_get_contents("$userDir/config/biopan_params.json");
$privateParams = json_decode($str, TRUE);
$srcFileName = $privateParams['srcFileName'];

# Load messages returned from parse_data.r script
$str = file_get_contents("$userDir/biopan/msg1.json");
$msg = json_decode($str, TRUE);
$na_values = $msg['na_values'][0];
$empty_rows = $msg['empty_rows'][0];
$error = $msg['error'][0];
$dup_lp = $msg['dup_lp'];
$str_dup_lp = join($dup_lp, ", ");


# Read summary.csv file to extract summary information
$summaryJson = file_get_contents("$userDir/biopan/summary.json");
$summaryParams = json_decode($summaryJson, TRUE);
# Total number of species
$totalSpecies = $summaryParams['total'];

# Species that are not recognized
$undefSpecies = $summaryParams['undef'];
# Species used for pathway analysis
$method = $summaryParams['pathway']['name'];
$processedSpecies = $summaryParams['pathway']['processed'];
$unprocessedSpecies = $summaryParams['pathway']['unprocessed'];
$sample_groups = $summaryParams['groups'];

$num_processed_species = sizeof($processedSpecies);
$num_unprocessed_species = sizeof($unprocessedSpecies);
$num_undef_species = sizeof($undefSpecies);
$num_sample_groups = sizeof($sample_groups);

# Get hostname URL of the server
$headers = apache_request_headers();
$hostURL = $headers["Host"];

include("includes_LM/top.php");
?>



<html>
<head>
    <title>LIPID MAPS® Lipidomics Gateway</title>
    <meta name="viewport" content="width=device-width, initial-scale=1  minimum-scale=1.0">
    <meta charset="utf-8">
    <meta http-equiv="x-ua-compatible" content="ie=edge">

    <link rel="stylesheet" href="css/style.css"/>
    <link rel="stylesheet" href="css/processing.css"/>
    <script type="text/javascript" src="js/librairies/jquery.min.js"></script>
    <script type="text/javascript" src="js/processing.js"></script>
</head>


<body>

    <?php
    # Error messages
    $message = "";

    if($error == true){
        $message = "LipidLynxX did not work properly. Check int the <a href='doc/step1#how-to-prepare-your-data.html' target='_blank'>documentation</a> that your file respects the format requested by BioPAN.";
    }else{
        if($sample_groups < 4){
            $message = "<li> The input file must contain at least 2 conditions and 2 samples per condition. </li>";
        }
        
        if($num_processed_species == 0 && $sample_groups > 4){
            $message = $message . "<li> There is no processed species. Look at the <a href='doc/step1#data-summary.html' target='_blank'>documentation</a> to understand why. </li>";
        }
    }

    if(!empty($message)){
        echo "  <div class='msg-box alert'>
                    <span class='closebtn' onclick='this.parentElement.style.display=\"none\";'>&times;</span>
                    <strong>Error:</strong>
                    <br>
                    <ul>
                        $message
                    </ul>
                </div>";
    }else{
        # Info messages
        $message = "";

        if($na_values == true){
            $message = "<li> There are NA value(s) in your input file. They have been replaced by 0. </li>";
        }

        if($empty_rows == true){
            $message = $message . " <li> There are empty row(s) in your input file. They have been deleted. </li>";
        }

        if(!empty($dup_lp)){
            $message = $message . "<li> There are duplicate lipids in your input file: $str_dup_lp. The values have been added. </li>";
        }

        if(!empty($message)){
            echo "  <div class='msg-box infos'>
                        <span class='closebtn' onclick='this.parentElement.style.display=\"none\";'>&times;</span>
                        <strong>Information:</strong>
                        <br>
                        <ul>
                            $message
                        </ul>
                    </div>";
        }
    }
    ?>

    <center><h1>BioPAN's Workflow</h1></center><br>


    <div id='anc_updates' style='line-height: 1.4;'>
        <?php
        echo "<br><strong>The data file <font color=#FF000>$srcFileName</font> has been successfully uploaded.</strong><br>";
        ?>
    </div>
    <br>
    <div id='summary'>
        <h4> Here is your data summary:</h4>
        <br>
        <table stype = 'margin-top:30px;' align="center">
            <tr>
                <th style="background-color: #b9c9fe; text-align: center;"></th>
                <th style="background-color: #b9c9fe; text-align: center;">Number</th>
            </tr>
            <tr>
                <td class="tooltips" style="background-color: #d2e4fc;">Unrecognised molecular species
                    <span class="tooltipstext">Species whose subclass is not recognised by BioPAN</span>
                </td>
                <td style="background-color: #e8edff; text-align: center;"><?php if ($num_undef_species > 0) {echo "<a target=\"_blank\" href='show_species.php?id=". $userID . "&type=undef'>" . $num_undef_species . "</a>";} else {echo $num_undef_species;} ?></td>
            </tr>
            <tr>
                <td class="tooltips" style="background-color: #d2e4fc;">Processed molecular species
                    <span class="tooltipstext">Species that are involved in at least one reaction</span>
                </td>
                <td style="background-color: #e8edff; text-align: center;"><?php if ($num_processed_species > 0) {echo "<a target=\"_blank\" href='show_species.php?id=". $userID . "&type=processed'>" . $num_processed_species . "</a>";} else {echo $num_processed_species;} ?></td>
            </tr>
            <tr>
                <td class="tooltips" style="background-color: #d2e4fc;">Unprocessed molecular species
                    <span class="tooltipstext">Species that are not involved in any reactions</span>
                </td>
                <td style="background-color: #e8edff; text-align: center;"><?php if ($num_unprocessed_species > 0) {echo "<a target=\"_blank\" href='show_species.php?id=". $userID . "&type=unprocessed'>" . $num_unprocessed_species . "</a>" ;} else {echo $num_unprocessed_species; }?></td>
            </tr>
            <tr>
                <td style="background-color: #d2e4fc;">Total</td>
                <td style="background-color: #e8edff; text-align: center;"><?php echo $totalSpecies;?></td>
            </tr>
        </table>
        <br>

        <tr>
            <td colspan=4 style='padding: 10px 0px 0px 0m'>There are guidelines on how interpret the data summary <a href='doc/step1.html#data-summary' target='_blank'>here</a>.</td>
        </tr>
    </div>
    <br>
    <br>
    <br>

    <?php
    if($num_processed_species > 0 && $sample_groups >= 4){
        echo "  <div id='anc_step2'>
                <h3>STEP 2: Assign conditions to different samples</h3>

                <div style='text-align: justify; text-justify: inter-word;'>In this step you have a chance to assign a condition to each sample.</div>
                <br>

                <form id='form' method='POST' onSubmit='submitForm();return false;'>
                <input type='hidden' id='userid' value='$userID'>
                <table id='assignTable' style='border-collapse: collapse;' align='center'>
                <thead>
                    <tr>
                        <th style='padding: 2px 10px 2px 10px; text-align: center; background-color: #b9c9fe;'>Sample name</th>
                        <th style='padding: 2px 0px 2px 10px; text-align: center; background-color: #b9c9fe; border-left: 2px solid white !important;'>Sample condition</th>
                    </tr> 
                </thead>
                <tbody>";

        $nb = 0;
        foreach ($sample_groups as $sample => $group) {
            echo "<tr>
                    <td style='padding:3px 10px 5px 0px; background-color:#e8edff;'>$sample:</td>
                    <input type='hidden' name='samples[]' value='$sample'>
                    <td id='columnTable' style='padding:3px 0px 5px 0px; text-align:center; background-color:#e8edff;'>
                        <input id='valuesTable$nb' type='text' name='groups[]' value='" . $group . "' style='width: 200px;'>
                    </td>
                </tr>";
            $nb += 1;
        }

        echo "      </tbody>
                <tr>
                <td colspan=2 style='padding:10px 0px 0px 0px; height:50px; text-align:center; vertical-align:bottom;'>";

        if($nb < 4){
            echo    "<input type='submit' disabled id='assignGroups' value='Assign conditions' style='cursor: no-drop; height: 30px; width: 150px; font-weight: bold;'>";
        }else{
            echo    "<input type='submit' id='assignGroups' value='Assign conditions' style='height: 30px; width: 150px; font-weight: bold; cursor: pointer;'>";
        }                
        
        echo "          </td>
                    </tr>
                </table>
                
                <script>
                // To control the cell color: red if less than 2 samples with the same condition name; otherwise green
                let valuesTable = document.querySelectorAll('[id^=\"valuesTable\"]');
                getColor();
                //Update when change the text box
                $('#assignTable td:last-child').on('change keydown paste input', function(){
                    getColor();
                });
                </script>

                <br>
                <br>

                <div class='overlay' style='display:none;'>
                    <div class='loading'>
                        <p>Processing data...</p>
                    </div>
                </div>

                <div class='less-than-one-sample-alert' style='display: none;'>
                    <span class='closebtn' onclick='this.parentElement.style.display='none';'>&times;</span> 
                    <center>Please introduce <strong>more than one sample condition.</strong></center>
                </div>
            
            </form>
        </div>";
    }
    ?>
</body>
</html>


<?php
include("includes_LM/bottom.php");
?>
