<?php
$userID = filter_input(INPUT_GET, 'id', FILTER_SANITIZE_STRING);
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
    $type = filter_input(INPUT_GET, 'type', FILTER_SANITIZE_STRING);
    # Read summary.csv file to extract summary information
    $summaryJson = file_get_contents("$userDir/biopan/summary.json");
    $summaryParams = json_decode($summaryJson, TRUE);
    $title = '';
    $data = array();
    
    # did LipidLynxX run?
    $lipidlynxx = $summaryParams['lipidlynxx'];

    # retrieve the classification of the lipid molecular species
    $undefSpecies = $summaryParams['undef'];
    $processedSpecies = $summaryParams['pathway']['processed'];
    $unprocessedSpecies = $summaryParams['pathway']['unprocessed'];

    $undefSpecies_dataset = $summaryParams['undef_dataset'];
    $processedSpecies_dataset = $summaryParams['processed_dataset'];
    $unprocessedSpecies_dataset = $summaryParams['unprocessed_dataset'];

    switch($type){
        case 'undef':
            $title = 'Unrecognised molecular species';
            $dataset = $undefSpecies_dataset;
            if($lipidlynxx == "no"){
                $dataset = $undefSpecies;
            }
            break;
        case 'unprocessed':
            $title = 'Unprocessed molecular species';
            $data = $unprocessedSpecies;
            $dataset = $unprocessedSpecies_dataset;
            if($lipidlynxx == "no"){
                $dataset = $unprocessedSpecies;
            }
            break;
        case 'processed':
            $title = 'Processed molecular species';
            $data = $processedSpecies;
            $dataset = $processedSpecies_dataset;
            if($lipidlynxx == "no"){
                $dataset = $processedSpecies;
            }
            break;
    }

include("includes_LM/top.php");
?>

<div id="show-species" style = "margin-top: 10px;">
    <h2><?php echo $title;?></h2>
    <br>
    <?php if($lipidlynxx == "yes" && $type != "undef"){
        echo "<p>Here is a table of the $title with a mapping between the nomenclature used by BioPAN and the nomenclature of the dataset. </p>";
    }
    ?>
    <table border="1" width="400">
        <thead>
            <tr>
                <?php 
                if($lipidlynxx == "yes" && $type != "undef"){
                    echo "<th>BioPAN nomenclature</th><th>Dataset nomenclature</th>";
                }else{
                    echo "<th>Lipid molecular species</th>";
                }
                ?>
            </tr>
        </thead>
        <tbody>
            <?php
            if($lipidlynxx == "yes" && $type != "undef"){
                if(count($dataset) == 1){
                    echo "<tr><td>" . $data[0] . "</td><td>" . $dataset . "</td></tr>";
                }
                else{
                    $x = 0;
                    foreach($dataset as $species){
                        echo "<tr><td>" . $data[$x] . "</td><td>" . $species . "</td></tr>";
                        $x += 1;
                    }
                }
            }else{
                if(count($dataset) == 1){
                    echo "<tr><td>" . $dataset . "</tr></td>";
                }
                else{
                    foreach($dataset as $species){
                        echo "<tr><td>" . $species . "</tr></td>";
                    }
                }
            }
            ?>
        </tbody>
    </table>
</div>
