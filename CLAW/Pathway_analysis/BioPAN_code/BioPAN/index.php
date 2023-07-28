<?php
$userID = substr(str_shuffle("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"), 0, 16);
# User's root folder
$userDir = "/lipidmaps/temp/$userID";
include("includes_LM/top.php");
?>

<script src='js/librairies/jquery.min.js'></script>
<link rel='stylesheet' href='https://fonts.googleapis.com/icon?family=Material+Icons'>

<center><h1>BioPAN: Bioinformatics Methodology For Pathway Analysis</h1></center><br>

<div style='text-align: justify; text-justify: inter-word; line-height: 1.2;'>

    <p>
        BioPAN for <strong>Bioinformatics Methodology For Pathway Analysis</strong> is a tool that allows users to upload their own mammalian lipidomics dataset 
        and perform a pathway analysis. 
        <br>
        In this analysis, you can explore systematic changes in lipid pathways at different levels: lipid subclass and lipid molecular species. 
        Those pathways will be highlighted and changes in gene activity will also be predicted.
    </p>
    <br>

    <p>BioPAN is fully integrated with LipidLynxX<sup><b>1</b></sup> to allow the user to upload dataset from different naming conventions and levels.</p>
    <button type="button" style='font-weight: bold; cursor:pointer; text-decoration:none;'>
                <a href="/lipidlynxx/" style="color: inherit; text-decoration: none" target="_blank">LipidLynxX</a>
            </button>
            - LipidLynxX is available on the LIPID MAPS website. <br/><br/>

    <br>
    <p>
        <span class="h3">Resources</span> <br/>
        <div class="pl-3">
            <button type="button" style='font-weight: bold; cursor:pointer; text-decoration:none;'>
                <a href="doc/index.html" style="color: inherit; text-decoration: none" target="_blank">Documentation</a>
            </button>
            - We recommend the user to read the detailed about BioPAN before using. <br/><br/>
            <button type="button" style='font-weight: bold; cursor:pointer; text-decoration:none;'>
                <a href="https://www.youtube.com/watch?v=EfTC_LuR3ak" style="color: inherit; text-decoration: none" target="_blank">Video Tutorial</a>
            </button>
            - Watch our tutorial video on YouTube.
        </div>
    </p>
    
    <br><br><br>
    <center>
        <button type="button" style='height: 35px; width: 160px; font-weight: bold; font-size: +1; text-decoration:none;'><a href="load_data.php" style="color: inherit; text-decoration: none">Start Analysis</a></button>
        <br><br><br><img src='doc/_images/workflow.png' height='430' border='1'>
    </center>

    <br><br>
    <p><sup><b>1</b></sup><strong> LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets.</strong> Zhixu Ni, Maria Fedorova. <br>
            <a href='https://www.biorxiv.org/content/10.1101/2020.04.09.033894v2' target='_blank'>bioRxiv, 2020.04.09.033894.</a></p>
        <br>

    <strong>Using lipidomics analysis to determine signalling and metabolic changes in cells.</strong>
    An Nguyen, Simon A Rudge, Qifeng Zhang, Michael JO Wakelam.</br>
    <a href='https://www.sciencedirect.com/science/article/pii/S0958166916302233' target='_blank'>Current Opinion in Biotechnology, Volume 43, 2017, Pages 96-103.</a></p>

</div>


<?php
include("includes_LM/bottom.php");
?>
