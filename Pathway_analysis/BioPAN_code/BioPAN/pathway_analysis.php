<?php
$userID = filter_input(INPUT_GET, 'id', FILTER_SANITIZE_STRING);
$userID = preg_replace("[^A-Za-z0-9]", "", $userID);
# User's root folder
$userDir = "/lipidmaps/temp/$userID";

# Show error if the page has been accessed with a wrong ID
if (empty($userID) || !file_exists($userDir)) {
    include("includes_LM/top.php");
    echo "<center><h1>BioPAN: Bioinformatics Methodology For Pathway Analysis</h1></center><br><br><br>
          <strong>ERROR:</strong> Unrecognised user \"$userID\".<br><br>
          Please go back to the <a href='index.php'>index</a> to get a new identifier.";
    include("includes_LM/bottom.php");
    exit;
}

# Load messages returned from process_data.r script
$str = file_get_contents("$userDir/biopan/msg2.json");
$msg = json_decode($str, TRUE);
$valid_groups = $msg['valid']['groups'];
$valid_freqs = $msg['valid']['freqs'];
$notvalid_groups = $msg['notvalid']['groups'];
$notvalid_freqs = $msg['notvalid']['freqs'];
$notvalid_str_group = join($notvalid_groups, ",");

$valid_lp_reaction = $msg['reaction']['lp'][0];
$valid_fa_reaction = $msg['reaction']['fa'][0];
$subset_pathway = $msg['subset']['pathway'][0];

$group_sample = $msg['group_sample']['equal'][0];
$is_dup = $msg['is_dup'][0];

//convert array into json
$str_group = join($valid_groups, ",");
$str_group_freq = join($valid_freqs, ",");
?>


<html>

<head>
    <title>LIPID MAPS® Lipidomics Gateway</title>
    <meta name="viewport" content="width=device-width, initial-scale=1  minimum-scale=1.0">
    <meta charset="utf-8">
    <meta http-equiv="x-ua-compatible" content="ie=edge">

    <link rel="icon" type="image/png" href="includes_LM/images/favicon.ico"/>
    <link rel="stylesheet" href="css/librairies/spliter.css">
    <link rel="stylesheet" href="css/librairies/tooltipster.bundle.min.css">   
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
    <link rel="stylesheet" href="css/librairies/jquery.cytoscape.js-panzoom.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <link rel="stylesheet" href="css/processing.css">
    <link rel="stylesheet" href="css/librairies/tippy.css">
    <link rel="stylesheet" href="css/librairies/bootstrap-duallistbox.css">
    <link rel="stylesheet" href="css/style.css">
    <link rel="stylesheet" href="css/index.css">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/jstree/3.2.1/themes/default/style.min.css" />
    
    <script src="js/librairies/jquery.min.js"></script>
    <script src="js/librairies/jquery-ui.js"></script>
    
    <script src="js/librairies/cytoscape.min.js"></script> <!-- version 3.10.2 -->
    <script src="https://unpkg.com/layout-base/layout-base.js"></script>
    <script src="js/librairies/cose-base.js"></script>
    <script src="js/librairies/cytoscape-cose-bilkent.js"></script>
    <script src="js/librairies/cytoscape-fcose.js"></script> 
    <script src="https://unpkg.com/popper.js"></script>
    <script src="js/librairies/cytoscape-popper.js"></script> <!-- version 1.0.2 -->
    <script src="js/librairies/jquery.cytoscape.js-panzoom.js"></script>
    
    <script src="js/librairies/FileSaver.js"></script>

    <script src="js/librairies/bootstrap.3.3.7.min.js"></script>
    <script src="js/librairies/jquery.bootstrap-duallistbox.js"></script>
    
    <script src="js/librairies/jstree.min.js"></script>
    <script src="js/librairies/splittedPanel.js"></script>
    <script src="js/librairies/tippy.all.js"></script>
    <script src="js/graph.js"></script>
</head>


<body onload="openAnalysisType(event, 'pathway', '<?php echo "$userDir"; ?>', '<?php echo "$str_group"; ?>','<?php echo "$str_group_freq"; ?>')">

    <!-- Popup window for graph exportation -->
    <div class="modal fade" id="modalSaveFile" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true" style="top:30%;">
        <div class="modal-dialog" role="document">
            <div class="modal-content">
            <div class="modal-header text-center">
                    <h5 class="modal-title">Export
                        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                            <span aria-hidden="true">&times;</span>
                        </button>
                    </h5>
                </div>
                <div class="modal-body mx-3">
                    <div class="md-form md-5">
                        <div class="btn-group btn-group-toggle" data-toggle="buttons">
                            <div class="custom-control custom-radio">
                                <input type="radio" class="custom-control-input" id="export-img" name="export-radios" value="graph" checked>
                                <label class="custom-control-label" for="export-img">Export graph</label>
                                <select id="export-graph">
                                    <option value="jpeg" selected>JPEG</option>
                                    <option value="png">PNG</option>
                                </select>
                            </div>
                            
                            <div class="custom-control custom-radio">
                                <input type="radio" class="custom-control-input" id="content-export" name="export-radios" value="content">
                                <label class="custom-control-label" for="content-export">Export graph content</label>
                                <select id="export-content">
                                    <option value="txt" selected>TXT</option>  
                                    <option value="json">JSON</option>              
                                </select>
                            </div>

                            <div class="custom-control custom-radio">
                                <input type="radio" class="custom-control-input" id="results-export" name="export-radios" value="results">
                                <label class="custom-control-label" for="results-export">Export results</label>
                                <select id="export-results">
                                    <option value="lipid-all" selected>All results</option>
                                    <option value="lipid-class">Lipid class results</option>
                                    <option value="lipid-species">Lipid species results</option> 
                                    <option value="fatty-acid">Fatty acid results</option>         
                                </select>
                            </div>
                        </div>  
                    </div>
                </div>
            
                <div class="modal-footer d-flex justify-content-center">
                    <button class="btn btn-default" id="export-btn">Export</button>
                </div>
            </div>
        </div>
    </div>
    <!--  end of popup window -->


    <!-- Popup window for graph legend -->
    <div class="modal fade" id="modalLegend" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true" style="top:15%;">
        <div class="modal-dialog" role="document">
            <div class="modal-content">
                <div class="modal-header text-center">
                    <h5 class="modal-title">Graph legend
                        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                            <span aria-hidden="true">&times;</span>
                        </button>
                    </h5>
                </div>
                <div class="modal-body mx-3 text-center">
                    <img id="myImg" src="doc/_images/Legend.png" alt="Graph legend" style="width:100%;">   
                </div>
            </div>
        </div>
    </div>
    <!--  end of popup window -->


    <!-- Popup window for edge click event -->
    <div class="modal fade" id="modalEdgeInfor" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true" style="top:20%;">
        <div class="modal-dialog" role="document">
            <div class="modal-content">
                <div class="modal-header text-center">
                    <h5 class="modal-title">Reaction table
                        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
                            <span aria-hidden="true">&times;</span>
                        </button>
                    </h5>
                </div>
                <div class="modal-body mx-3" id="table-infor">
                    <div class="modal-header text-center">
                        <h5 class="modal-title">Reactant group</h5>
                        <select multiple="multiple" size="10" name="duallistbox_reactant" title="duallistbox_reactant" id="dual-list-box-re">
                        </select>
                    </div>
                    <div class="modal-header text-center">
                        <h5 class="modal-title">Product group</h5>
                        <select multiple="multiple" size="10" name="duallistbox_product" title="duallistbox_product" id="dual-list-box-pro">
                        </select>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <!--  end of popup window -->


    <div id="box">
        <a id="leftbox" href="/"><img src="/includes_LM/images/logo_R_transparent.png" width="60" alt="logo"><span style="font-size:2.1rem; padding-left:0.8rem;">LIPID MAPS® Lipidomics Gateway</span></a>
        <p id="middlebox">BioPAN: pathways analysis</p>
    </div>


    <?php
    # Error messages
    $message = "";

    if(count($valid_groups) < 2){
        $message = "<li> The input file must contain at least 2 conditions and each group must contain at least 2 samples. </li>";
    }

    if($is_dup == false){
        $message = "<li> The input file must contain at least 2 conditions with the same number of samples. </li>";
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
        
        if(count($notvalid_groups) > 0){
            $message = "<li> The conditions(s) '$notvalid_str_group' have been removed. Each condition must contain at least 2 samples. </li>";
        }

        if($valid_lp_reaction == false && $valid_fa_reaction == false){
            $message = $message . "<li> No reactions found. </li>";
        }else{
            if($valid_lp_reaction == false ){
                $message = $message . "<li> No lipid reactions found. </li>";
            }
            if($subset_pathway == false){
                $message = $message . "<li> No reactions found for the pathway subset of lipid data. </li>";
            }
            if($valid_fa_reaction == false ){
                $message = $message . "<li> No fatty acid reactions found. </li>";
            }
        }

        if($group_sample == false){
            $message = $message . "<li> Not all conditions have the same number of samples. The conditions to be compared must have the same number of samples. </li>";
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


    <div style="float:right; margin-right:0.6%;">
        <button id="btnSave" type="button" title="Export" data-toggle="modal" data-target="#modalSaveFile">
            Export
        </button>
        <button type="submit" title="Documentation" onclick="window.open('doc/step3_build.html','_blank')">
            Documentation
        </button>
    </div>
  
   

    <div id="pathway" class="tabcontent" style="overflow=auto;">
        <div class="pretty-split-pane-frame">
            <div class="split-pane fixed-left" id="split-pane-1">
                <div class="split-pane-component" id="left-component">
                  
                    <!-- Graph legend panel -->
                    <div id="graph-panel" class="panel panel-info  img-rounded panel-hide">
                        <div class="panel-heading">Graph legend</div>            
                        <div class="panel-body">
                            <button type="button" title="Legend" data-toggle="modal" data-target="#modalLegend">
                                <img id="myImg" src="doc/_images/Legend.png" alt="Graph legend" style="width:100%;">   
                            </button>
                        </div>
                    </div>
                   

                    <!-- Pathway options panel -->
                    <div id="display-panel" class="panel panel-info img-rounded panel-hide">
                        <div class="panel-heading">Pathway options</div>            
                        <div class="panel-body">
                            <label for="groups" class="tooltip u-pull-left" title="Choose conditions for comparison">Group comparison</label>
                            <div class="u-cf"></div>
                            <div>
                                Groups for comparison:
                                <div class="row" style="line-height:10px;">
                                    <div class="col-sm-6">
                                        <div class="form-group">
                                            <label for="disease_group">Condition of interest</label>
                                            <select id="disease_group" class="form-control" style="height:auto;">
                                                <option value="group1" selected>Group 1</option>
                                                <option value="group2">Group 2</option>
                                                <option value="group3">Group 3</option>
                                            </select>
                                        </div>
                                    </div>
                                    <div class="col-sm-6">
                                        <div class="form-group">
                                            <label for="control_group">Control condition</label>
                                            <select id="control_group" class="form-control" style="height:auto;">
                                                <option value="group1" selected>Group 1</option>
                                                <option value="group2">Group 2</option>
                                                <option value="group3">Group 3</option>
                                            </select>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            <label for="pathway_type" class="tooltip u-pull-left" title="Switch to a different pathway">Pathway type</label>
                            <div class="u-cf"></div>
                            Type: 
                            <select id="pathway_type">
                                <option value="lp" selected>Lipid</option>
                                <option value="fa">Fatty acid</option>
                            </select>
                            
                            <label for="pathway_status" class="tooltip u-pull-left" title="Switch to pathway status">Pathway status</label>
                            <div class="u-cf"></div>
                            Status: 
                            <select id="pathway_status">
                                <option value="active" selected>Active</option>
                                <option value="most_active" >Most active</option>
                                <option value="suppressed">Suppressed</option>
                                <option value="most_suppressed">Most suppressed</option>
                            </select>

                            <div id="pw_level">
                                <label for="pathway_level" class="tooltip u-pull-left" title="Switch to pathway level">Pathway level</label>
                                <div class="u-cf"></div>
                                Level: 
                                <select id="pathway_level">
                                    <option value="class" selected>Lipid subclass</option>
                                    <option value="species">Lipid molecular species</option>
                                </select>
                            </div> 
                            
                            <div id = "filter" style="margin-top:10px;">
                                Filter: <input type="text" id="search" placeholder="search"><button id="searchBox">Search</button>
                                <div id="class_reaction" style="display:none;">
                                </div>
                                <div id="species_reaction" style="display:none;">
                                </div>
                                <div id="class_pathway" style="display:none;">
                                </div>
                                <div id="species_pathway" style="display:none;">
                                </div>
                            </div>

                            <div id="sub_set" style="margin-top: 15px;">
                                <label for="subset" class="tooltip u-pull-left" title="Select subset of data">Select subset of data</label>
                                <div class="u-cf"></div>
                                Subset of lipid data: 
                                <select id="subset">
                                    <option value="reaction" selected>Reactions</option>
                                    <option value="pathway" >Pathways</option>
                                </select>
                            </div>
                        </div>
                    </div>

                    
                    <!-- Pathway calculation panel -->
                    <div id="display-panel" class="panel panel-info img-rounded panel-hide">
                        <div class="panel-heading">Pathway calculation</div>
                            <div class="panel-body">
                                <label for="pvalue" class="tooltip u-pull-left" title="Choose a p-value">P-value</label>
                                <div class="u-cf"></div>
                                P-value: 
                                <select id="pvalue">
                                    <option value="0.05" selected>0.05</option>
                                    <option value="0.02">0.02</option>
                                    <option value="0.01">0.01</option>
                                </select>
                                
                                <label for="type_of_expr" class="tooltip u-pull-left" title="Choose type of experimental measurement">Paired data</label>
                                <div class="u-cf"></div>
                                Paired data: 
                                <select id="type_of_expr">
                                    <option value="FALSE" selected>No</option>
                                    <option value="TRUE">Yes</option>
                                </select>
                            </div>
                            <div style="text-align: center;">
                                <input type="button" id="calculate" value="Calculate pathways" style="height:30px; width:150px; font-weight:bold; margin-top:5px; margin-bottom:10px;">
                            </div>
                        </div>
                    </div>
                        
                
                    <!--verical divider-->
                    <div class="split-pane-divider" id="vertical-divider"></div>
                    <!--verical divider-->

                    <div class="split-pane-component" id="right-component">
                        <div class="split-pane fixed-bottom" id="split-pane-2">

                            <div class="split-pane-component" id="top-component-2">                               
                                <div class="pretty-split-pane-component-inner">
                                
                                    <div class="diagram-section component-section" style="height: 350px;">
                                        <div id="cy-container">
                                            <div id="cy"></div>
                                            
                                            <div class="graph-loading" style="display:none;">
                                                <span class="fa fa-refresh fa-spin" style="font-size: 50px;"></span>
                                            </div>

                                        </div>
                                        
                                    </div>
                                </div>
                            </div>
                            
                            <!-- Horizontal splitter -->
                            <div class="split-pane-divider" id="horizontal-divider-2"></div>
                            <!-- Horizontal splitter -->
                            
                            <!-- Table panel -->
                            <div class="split-pane-component" id="bottom-component-2">
                                <div class="tab">
                                    <button class="tabresults" id="all-lp" style="text-decoration: none" onclick="openAnalysisType1(event, 'lipid-all', '<?php echo "$userDir"; ?>')">All results
                                    <button class="tabresults" id="lp-cl" style="text-decoration: none;" onclick="openAnalysisType1(event, 'lipid-class', '<?php echo "$userDir"; ?>')">Lipid class results</button>
                                    <button class="tabresults" id="lp-sp" style="text-decoration: none;" onclick="openAnalysisType1(event, 'lipid-species', '<?php echo "$userDir"; ?>')">Lipid species results</button>
                                    <button class="tabresults" id="all-fa" style="text-decoration: none;" onclick="openAnalysisType1(event, 'fatty_acid', '<?php echo "$userDir"; ?>')">Fatty acid results</button>

                                </div>
                                    
                                <div id="lipid-all" class="tabdetail inactive">
                                    <div class="split-pane-bottom-component" id="bottom-component-3">
                                        <div class="pretty-split-pane-component-inner" id="lp-all">
                                        </div>
                                    </div>
                                </div>
                                
                                <div id="lipid-class" class="tabdetail inactive">
                                    <div class="split-pane-bottom-component" id="bottom-component-3">
                                        <div class="pretty-split-pane-component-inner" id="lp-class">
                                        </div>
                                    </div>
                                </div>

                                <div id="lipid-species" class="tabdetail inactive">
                                    <div class="split-pane-bottom-component" id="bottom-component-3">
                                        <div class="pretty-split-pane-component-inner" id="lp-species">
                                        </div>
                                    </div>
                                </div>
                                
                                <div id="fatty_acid" class="tabdetail inactive">
                                    <div class="split-pane-bottom-component" id="bottom-component-3">
                                        <div class="pretty-split-pane-component-inner" id="fa">
                                        </div>
                                    </div>
                                </div>

                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</body>
</html>