window.is_root_selected = 1;
window.userDir;
window.is_refresh = 0;
window.tippyNode = null;
window.all_re_arr = [];
window.all_pro_arr = [];
window.listBoxRe = null;
window.listBoxPro = null;
window.cy = null;


//Handle refresh event
if(performance.navigation.type == 1){
    is_refresh = 1;
}


function openAnalysisType1(event, name, userDir){
    var event = event || window.event;
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName('tabdetail');
    for(i = 0; i < tabcontent.length; i++){
        tabcontent[i].style.display = 'none';
    }
    tablinks = document.getElementsByClassName('tabresults');
    for(i = 0; i < tablinks.length; i++){
        tablinks[i].className = tablinks[i].className.replace('active', '');
    }
    document.getElementById(name).style.display = 'block';
    event.target.className += ' active';
}


function fillGroupOptions(selectorId, groupArr, selectedValue){
    let dropdown = document.getElementById(selectorId);
    // optional: clear all existing options first:
    if(dropdown){
        dropdown.innerHTML = "";
        // populate list with options:
        for(let i = 0; i < groupArr.length; i++) {
            let opt = groupArr[i];
            if(opt == selectedValue){
                dropdown.innerHTML += '<option value=\'' + opt + '\' selected>' + opt + '</option>';    
            }else{
                dropdown.innerHTML += '<option value=\'' + opt + '\'>' + opt + '</option>';    
            }
        }
    }
}



// **** Functions related to the export of the results tables **** //

//Get status to write it on the title and file title for export
function status(strPwStatus){
    let title, fileTitle;
    switch(strPwStatus){
        case 'active':
            title = fileTitle = 'Active';
            break;
        case 'most_active':
            title = 'The most active pathways';
            fileTitle = 'MostActive';
            break;
        case 'suppressed':
            title = fileTitle = 'Suppressed';
            break;
        case 'most_suppressed':
            title = 'The most suppressed';
            fileTitle = 'MostSuppressed';
    }
    return [title, fileTitle];
}


//Get level and status from a json file and reformate the information
function levelStatus(json_path){
    let levels = '';
    let status = '';
    let file = json_path.split('/');
    file = file[file.length-1].split('_');

    if(file[0] == 'lp'){
        levels = 'lipid ' + file[1];
    }else{
        levels = 'fa';
    }

    if(json_path.includes('most_suppressed')){
        status = 'most suppressed';
    }else if(json_path.includes('most_active')){
        status = 'most active';
    }else if(json_path.includes('suppressed')){
        status = 'suppressed';
    }else{
        status = 'active';
    }

    return [levels,status]
}

//Fill the results table to export (text format)
function fillFile(res, subset, type){
    let genes;
    let href = '';
    let write_genes = '';
    genes = res['data']['gene'].split(',');
    fLen = genes.length;
    if (genes != ''){
        if(fLen == 1 && genes[0] == 'NA'){ // no genes associated to a reaction
            href = 'No genes have yet been identified'
        }
        for(let j = 0; j < fLen; j++) {
            if(genes[j].trim() !== 'NA'){
                write_genes += genes[j] + ', '
                href += 'https://lipidmaps.org/data/proteome/LMPD_table.php?GENE_SYMBOL=' + genes[j].trim() + ' ; ';
            }
        }
        write_genes = write_genes.substring(0, write_genes.length-2);
        href = href.substring(0, href.length-2);
    }

    let pathway_list = res['data']['pathway'].split(/[;&]/);
    let pathway = '';
    for(let j in pathway_list) {
        if (j%2 == 0){
            pathway += pathway_list[j] + '->';
        }
    }
    pathway = pathway.substring(0, pathway.length-2);
    
    output = pathway + '\t';
    if(subset == 'pathway'  && type != 'fa'){
        output += res['data']['class'] + '\t';
    }
    output += res['data']['score'] + '\t' + write_genes + '\t' + href + '\n';
    return output;
}


//Save the requested results table
function saveTables(output, fileTitle, name, subset, diseaseGr, controlGr, pValue){
    let fileName = subset + '_' + diseaseGr + '_vs_' + controlGr + '_' + pValue + '.tsv';
    if(name == 'lipid-all'){
        fileName = 'All_' + fileName;
    }else if(name == 'lipid-subclass'){
        fileName = 'LipidSubclass_' + fileTitle + '_' + fileName;
    }
    else if(name == 'lipid-species'){
        fileName = 'LipidMolecularSpecies_' + fileTitle + '_' + fileName;
    }
    else{
        fileName = 'FattyAcid_' + fileTitle + '_' + fileName;
    }
    fileToSave = new File([output], {
        type: 'text/plain',
        name: fileName
    });
    saveAs(fileToSave, fileName);
}


//Export a results table
function exportResultTables(json_path, name, subset, strPwStatus, diseaseGr, controlGr, pValue){
    return $.ajax({
        type: 'GET',
        url: 'get_json_data_table.php', 
        data: {path: json_path},

        success: function(result){
            let fileTitle = status(strPwStatus)[1];
            let output = '';

            if(result[0] == 'null' && result.length == 1){
                output = name.replace('-', ' ').replace('_', ' ') + ' ' + status(strPwStatus)[0].toLowerCase();
                output += ' ' + subset + ' (' + diseaseGr + ' vs ' + controlGr + ')' + '\n';
                output += 'No results are found.' + '\n';
            
            }else{
                for(let k = 0; k < result.length; k++){
                    let type = '';
                    let levelsStatusP = levelStatus(json_path[k]);
                    let levels = levelsStatusP[0];
                    let statusP = levelsStatusP[1];

                    if(result.length == 1){
                        output = name.replace('-', ' ').replace('_', ' ') + ' ' + status(strPwStatus)[0].toLowerCase();
                        type = name;
                    }else{
                        output += levels + ' ' + statusP;
                        type = levels;
                    }

                    output += ' ' + subset + ' (' + diseaseGr + ' vs ' + controlGr + ')' + '\n';

                    if(result[k] == 'null'){
                        output += 'No results are found.' + '\n' + '\n';
                    
                    }else{
                         if(subset == 'pathway' && type != 'fa'){
                            output += 'Pathways';
                            output +=  '\t' + 'Classification(s) of the pathways';
                        }else{
                            output += 'Reactions chains';
                        } 

                        output += '\t' + 'Z-score' + '\t' + 'Predicted genes' + '\t' + 'Predicted genes links' + '\n';
                        res = JSON.parse(result[k])['pathways'];
                        
                        for(let i in res) {
                            output += fillFile(res[i], subset, type);
                        }
                        output += '\n';
                    }
                }
            }
            saveTables(output, fileTitle, name, subset, diseaseGr, controlGr, pValue);
        },
        error: function(xhr, status, err) { alert('Ajax request failed'); }
    });
}



// **** Functions related to the results tables **** //

//Get the correct title for the results tables and the type of the reaction
function getInfo(json_path, result, name, fileTitle, subset, diseaseGr, controlGr){
    let levelsStatusP = levelStatus(json_path);
    let levels = levelsStatusP[0];
    let statusP = levelsStatusP[1];
    if(result.length == 1){
        title = name.replace('-', ' ').replace('_', ' ') + ' ' + fileTitle;
        type = name;
    }else{
        title = levels.replace('class','subclass') + ' ' + statusP;
        type = levels;
    }
    
    title += ' ' + subset + ' (' + diseaseGr + ' vs ' + controlGr + ')';
    return [title, type];
}


//Write the title if there are no results
function noRes(title){
    output = '<div class="title-section">';
    output += '<p>' + title + '</p>' + '</div>';
    output += '<div class="section-content">';
    output += '<p style="font-size:16px;">No results are found.</p>';
    output += '</div>';
    return output
}


//Write the title of the results tables and the columns names
function headTable(title, subset, type){
    let output = '<div class="title-section">';
    output += '<p>' + title + '</p>' + '</div>';
    output += '<div class="section-content">';

    output += '<table class="table table-striped table-bordered" style="width:100%;"';
    output += '<thead><tr>';
    
    if(subset == 'pathway' && type != 'fa'){
        output += '<th>Pathways</th>';
        output += '<th>Classification(s) of the pathways</th>';
    }else{
        output += '<th>Reactions chains</th>';
    }

    output += '<th>Z-score</th>';
    output += '<th>Predicted genes</th>';
    output += '</tr><thead><tbody>';

    return output;
}


//Fill the results tables (HTML format)
function fillTable(res, title, subset, type){
    let genes;
    let href = '';
    genes = res['data']['gene'].split(',');
    fLen = genes.length;

    if(fLen == 1 && genes[0] == 'NA'){ // no genes associated to a reaction
        href = 'No genes have yet been identified'
    }else{
        for(let j = 0; j < fLen; j++) {
            if(genes[j].trim() !== 'NA'){
                href += "<a href = 'https://lipidmaps.org/data/proteome/LMPD_table.php?GENE_SYMBOL=" + genes[j].trim() + "' target='_blank'> <i>" + genes[j] + "</i></a>, ";
            }
        }
        href = href.substring(0, href.length-2);
    }
    
    let path = '';
                                
    if(title.includes('lipid species') && title.includes('reaction')){
        let paths = res['data']['pathway'].split('&#8594;');
        
        for(let j = 0; j < paths.length - 1; j++){
            path += paths[j] + ' &#8594; ';
        }
        path = path + paths[paths.length - 1];
    }else{
        path = res['data']['pathway'];
    }
    
    let output = '<tr><td>' + path + '</td><td>';
    
    if(subset == 'pathway' && type != 'fa'){
        output += res['data']['class'] + '</td><td>';
    }

    output += res['data']['score'] + '</td><td>' + href + '</td></tr>';
    return output;
}


//Visualise the requested results table
function visualiseTables(name, output){
    if(name == 'lipid-all'){
        $('#lp-all').html(output);
    }else if(name == 'lipid-subclass'){
        $('#lp-class').html(output);
    }else if(name == 'lipid-species'){
        $('#lp-species').html(output);
    }else{
        $('#fa').html(output);
    }
}


//Update the results tables: formatting the json result file to HTML
function updateResultTables(json_path, name, subset, pwStatus, diseaseGr, controlGr){
    return $.ajax({
        type: 'GET',
        url: 'get_json_data_table.php', 
        data: {path: json_path},

        success: function(result){
            let output = '';
            let title = '';

            let fileTitle = status(pwStatus)[0].toLowerCase();

            if(result[0] == 'null' && result.length == 1){
                title = name.replace('-', ' ').replace('_', ' ') + ' ' + fileTitle;
                title += ' ' + subset + ' (' + diseaseGr + ' vs ' + controlGr + ')' ;
                output = noRes(title);
            }else{
                for(let k = 0; k < result.length; k++){
                    let infos = getInfo(json_path[k], result, name, fileTitle, subset, diseaseGr, controlGr);
                    title = infos[0];
                    let type = infos[1];
                    
                    if(result[k] == 'null'){
                        output += noRes(title);
                    }else{
                        output += headTable(title, subset, type);
                        res = JSON.parse(result[k])['pathways'];
                        for(let i = 0; i < res.length; i++) {
                            output += fillTable(res[i], title, subset, type);
                        }
                        output += '</tbody></table>';
                    }
                }
            }
            output += '</div>';

            visualiseTables(name, output);
        },
        error: function(xhr, status, err) { alert('Ajax request failed'); }
    });
}



//Fill the data tree
function fillDataTree(){
    let pwLevel = ['class','species'];
    let subset = ['reaction','pathway'];
    for(var i = 0; i < 2; i++){
        for(var j = 0; j < 2; j++){  
            let dataPath = userDir + '/biopan/lp_' + pwLevel[i] + '_' + subset[j] + '.json';
            let dataTree = '#' + pwLevel[i] + '_' + subset[j];
            $.ajax({
                type: "GET",
                url: "get_json_data.php", 
                data: {path: dataPath},
                dataType: "json",
                success: function(res){
                    if(res != null){
                        let data = [{
                            'text' : 'Select all',
                            'id' : 'root',
                            'state' : {'opened': false, 'selected': true},
                            children: res
                        }];
                        $(dataTree).jstree({
                            'plugins': ['search', 'checkbox', 'wholerow'],
                            'core': {
                                "check_callback": true,
                                'data': data,
                                'animation': false,
                                'themes': {'icons': false}
                            },
                            'search': {'show_only_matches': true, 'show_only_matches_children': false},
                            'checkbox': {'two_state': false, 'whole_node':false, },
                        }).on('ready.jstree', function (e, data, selected, event) {
                            var ref = $(dataTree).jstree(true);
                            ref.close_all();
                        });
                    }
                }
            });
        }
    }
}


//Get the title for the "isPaired" option (TRUE: paired and FALSE:notpaired)
function getPairedTitle(isPaired){
    let paired = 'paired';
    if(isPaired == 'FALSE'){
        paired = 'notpaired';
    }
    return paired;
}


//Get the prefix path of the graph to display
function getGraphDatasetPrefix(){
    let subset = document.querySelector('#subset').value;
    let pwType = document.getElementById('pathway_type').value;
    let pwStatus = document.getElementById('pathway_status').value;
    let diseaseGr = document.getElementById('disease_group').value;
    let controlGr = document.getElementById('control_group').value;
    let isPaired = document.getElementById('type_of_expr').value;
    let paired = getPairedTitle(isPaired);
    
    let pathway_level='';
    let prefix_arr;

    if(pwStatus.indexOf('suppressed') >= 0){
        pwStatus = 'suppressed';
    }else{
        pwStatus = 'active';
    }
    if(pwType == 'lp'){
        pathway_level= document.getElementById('pathway_level').value;
        prefix_arr = ['lp', pathway_level, subset, diseaseGr, controlGr, pwStatus, paired];
    }else{
        prefix_arr = ['fa', diseaseGr, controlGr, pwStatus, paired];
    }
    let prefix_str = prefix_arr.join("_");

    return prefix_str;
}


//Get the full path of the graph to display
function getGraphDatasetPath(){
    let json_file_name = getGraphDatasetPrefix().concat(".json");
    let json_path = [userDir, 'biopan', json_file_name].join("/");
    return json_path;
}

//Get the full path of the filter graph to display
function getGraphDatasetPathFilter(){
    let json_file_name = getGraphDatasetPrefix().concat("_filter.json");
    let json_path = [userDir, 'biopan', json_file_name].join("/");
    return json_path;
}


//Get the path of the results tables
function getSigPathwayTablePath(name){
    let subset = document.querySelector('#subset').value;
    let pwStatus = document.getElementById('pathway_status').value;
    let diseaseGr = document.getElementById('disease_group').value;
    let controlGr = document.getElementById('control_group').value;
    let pValue = document.querySelector("#pvalue").value;
    let isPaired = document.getElementById('type_of_expr').value;
    let paired = getPairedTitle(isPaired);
    
    if(name == 'lipid-all'){
        let json_path = [userDir, 'biopan'].join("/");
        levels = ["lp_class","lp_species","fa"];
        statusP = ["active","most_active","suppressed","most_suppressed"];
        let lp_all_json_path = [];
        for(let i = 0; i < levels.length; i++){
            for(let j = 0; j < statusP.length; j++){
                let path = json_path + '/' + levels[i] + '_';
                if(levels[i] != "fa"){
                    path += subset + '_';
                }
                path += diseaseGr + '_' + controlGr + '_' + statusP[j] + '_' + pValue + '_' + paired + '_tbl.json';
                lp_all_json_path.push(path); 
            }
        }
        return lp_all_json_path;

    }else if(name == 'lipid-subclass'){
        prefix_arr = ['lp_class', subset, diseaseGr, controlGr, pwStatus, pValue, paired, "tbl"];
    }else if(name == 'lipid-species'){
        prefix_arr = ['lp_species', subset, diseaseGr, controlGr, pwStatus, pValue, paired, "tbl"];    
        
    }else{
        prefix_arr = ['fa', diseaseGr, controlGr, pwStatus, pValue, paired, "tbl"];
    }

    let prefix_str = prefix_arr.join("_");
    let json_file_name = prefix_str.concat(".json");
    let json_path = [userDir, 'biopan', json_file_name].join("/");
    return [json_path];
}


//Get the path containing the nodes and edges to highlight (= with status)
function getHighLightedGraphPath(){
    let pwType = document.getElementById('pathway_type').value;
    let pwStatus = document.getElementById('pathway_status').value;
    let subset = document.querySelector('#subset').value;
    let diseaseGr = document.getElementById('disease_group').value;
    let controlGr = document.getElementById('control_group').value;
    let pValue = document.getElementById('pvalue').value;
    let isPaired = document.getElementById('type_of_expr').value;
    let paired = getPairedTitle(isPaired);

    let pathway_level;
    let prefix_arr;
    if(pwType == 'lp'){
        pathway_level= document.getElementById('pathway_level').value;
        prefix_arr = ['lp', pathway_level, subset, diseaseGr, controlGr, pwStatus, pValue, paired];
    }else{
        prefix_arr = ['fa', diseaseGr, controlGr, pwStatus, pValue, paired];
    }

    let prefix_str = prefix_arr.join("_");
    let json_file_name = prefix_str.concat(".json");
    let json_path = [userDir, 'biopan', json_file_name].join("/");
    
    return json_path;
}


// Get the type (graph or table) and result of file to save
function getFileTypeToSave(){
    let radios = document.getElementsByTagName('input');
    let radioValue;
    let res;
    let type = 'graph';
    for(let i = 0; i < radios.length; i++) {
        if(radios[i].type === 'radio' && radios[i].checked) {
            // get value, set checked flag or do whatever you need to
            radioValue = radios[i].value;
        }
    }   
            
    switch(radioValue){
        case 'graph':
            res = document.getElementById('export-graph').value;
            break;
        case 'content':
            res = document.getElementById('export-content').value;
            break;
        case 'results':
            res = document.getElementById('export-results').value;
            type = 'table';
    }
    return [res,type];
}


//Update the 4 results tables to match with the requested results
function updateTables(subset, pwStatus, diseaseGr, controlGr){
    let lp_all_json_path = getSigPathwayTablePath('lipid-all');
    let lp_class_json_path = getSigPathwayTablePath('lipid-subclass');
    let lp_species_json_path = getSigPathwayTablePath('lipid-species');
    let fa_json_path = getSigPathwayTablePath('fa');

    let c1 = updateResultTables(lp_all_json_path, 'lipid-all', subset, pwStatus, diseaseGr, controlGr);
    let c2 = updateResultTables(lp_class_json_path, 'lipid-subclass', subset, pwStatus, diseaseGr, controlGr);
    let c3 = updateResultTables(lp_species_json_path, 'lipid-species', subset, pwStatus, diseaseGr, controlGr);
    let c4 = updateResultTables(fa_json_path, 'fa', subset, pwStatus, diseaseGr, controlGr);
    
    $.when(c1, c2, c3, c4).then(function (a1, a2, a3, a4) {
    }, function (jqXHR, textStatus, errorThrown) {
        let x1 = c1;
        let x2 = c2;
        let x3 = c3;
        let x4 = c4;
        if (x1.readyState != 4){
            x1.abort();
        }
        if (x2.readyState != 4){
            x2.abort();
        }
        if (x3.readyState != 4){
            x3.abort();
        }
        if (x4.readyState != 4){
            x4.abort();
        }
        alert('Either c1, c2, c3 or c4 failed!');
    });
}


// Final functions
function openAnalysisType(event, name, userDir, str_groups, str_group_freq) {
    let i, tabcontent, tablinks;
    
    tabcontent = document.getElementsByClassName('tabcontent');
    for(i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = 'none';
    }
    tablinks = document.getElementsByClassName('tablinks');
    for(i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace('active', '');
    }
    document.getElementById(name).style.display = 'block';
    event.target.className += ' active';
    
    this.userDir = userDir;
    let elem = document.getElementById('cy');
    // fill group options
    let groups = str_groups.split(',');
    let freqs = str_group_freq.split(',');
    let gr_freq = new Array();
    for(let i = 0; i < groups.length; i++){
        gr_freq[groups[i]] = freqs[i];
    }

    fillGroupOptions('control_group', groups, groups[0]);
    fillGroupOptions('disease_group', groups, groups[1]);

    let pwLevel = document.getElementById('pathway_level');
    let pwType = document.getElementById('pathway_type');
    let subset = document.getElementById('subset');
    let pwStatus = document.getElementById('pathway_status');
    let pValue = document.getElementById('pvalue');
    let diseaseGr = document.getElementById('disease_group');
    let controlGr = document.getElementById('control_group');
    let isPaired = document.getElementById('type_of_expr');

    let isCondition = false; // disease & control
    let isOther = false; // other: pValue and isPaired

    // fill data tree
    fillDataTree();

    // display data tree
    let dtKey = '#' + pwLevel.value + '_' + subset.value;


    //Filter the pathways with the selected lipid
    function filterPathway(selectedItemStr){
        let path = getGraphDatasetPath();
        $.ajax({
            type: 'GET',
            url: 'filter_pathway.php', 
            data: {path:path, selectedItems:selectedItemStr, userDir:userDir},
            dataType: 'text',
            success: function(){
                $('.graph-loading').hide();
                if(selectedItemStr.length > 0){
                    Promise.resolve(userDir).then(getGraphDatasetPathFilter).then(applyGraph);
                    updateTables(subset.value, pwStatus.value, diseaseGr.value, controlGr.value);
                }else{
                    cy.elements().remove();
                }
            }
        });
    }


    //On new search
    function onSearch(value){
        $('.graph-loading').show();
        if(to) {clearTimeout(to);}
        to = setTimeout(function(){
            let toRemove = ["Select all", "Glycerolipids and Glycerophospholipids", "Sphingolipids", "Triglycerol metabolism"];
                        
            let selectedItems = $(dtKey).jstree("search", value)[0].childNodes[0].innerText;
            selectedItems = selectedItems.split(/\r?\n/);
            selectedItems = selectedItems.filter((el) => !toRemove.includes(el));
            selectedItems = selectedItems.filter((el) => /\S/.test(el));
            let selectedItemStr = selectedItems.join(",");
            if(selectedItemStr.length < 1){
                drawGraph();
                $(dtKey).jstree("check_all");
                 $('.graph-loading').hide();
            }else{
                filterPathway(selectedItemStr);
            }            
        }, 250);
    }


    //Filter lipids after clicking on the search button
    let searchBtn = document.getElementById('searchBox');
    searchBtn.addEventListener('click', function(){
        let value = $("#search").val();
        onSearch(value);
    });

    
    //Get the name of the selected lipids in the tree
    function getSelectedNameTree(){
        let data = $(dtKey).jstree('get_selected', true);
        let lp = [];
        for(i = 0; i < data.length; i++){
            lp.push(data[i].text);
        }
        return lp;
    }

    //Get the name of the lipids from ids of the tree
    function getCorrespondanceTree(child){
        let data = $(dtKey).jstree('get_selected', true);
        lp = []
        for(i = 0; i < data.length; i++){
            if(child.indexOf(data[i].id) > -1){
                lp.push(data[i].text);
            }
        }
        return lp;
    }

    //Update filter tree when there is already a filter
    function filterOnSearch(data){
        let tree_path = [userDir, "biopan", "selected_items.csv"].join('/');
        $.ajax({
            type: 'GET',
            url: 'get_csv_data.php', 
            data: {path: tree_path},
            dataType: 'text',
            success: function(result){
                res = result.split(",");
                res[res.length - 1] = res[res.length - 1].replace(/(\r\n|\n|\r)/gm,"");
                let bool = false;
                if(data.action == "deselect_node"){ // if deselected node: remove it from the results
                    bool = true;
                    let child = data.node.children;
                    if(child.length > 0){ // if the deselected node(s) are not the last children
                        let match = getSelectedNameTree();
                        res = res.filter(value => match.includes(value));
                    }
                    else{
                        let index = res.indexOf(data.node.text);
                        if (index > -1) {
                            res.splice(index, 1);
                        }
                    }
                }
                if(data.action == "select_node"){ // if new selected node: add it to the results
                    bool = true;
                    let child = data.node.children;

                    if(child.length > 0){ // if the selected node(s) are not the last children
                        let newLp = getCorrespondanceTree(child);
                        res.push(newLp);
                    }else{
                        res.push(data.node.text);
                    }
                }
                if(bool){
                    let selectedItemStr = res.join(",");
                    filterPathway(selectedItemStr);
                }
            }
        });
    }


    //Check if it is a new search
    function defineFilter(data){
        let objects = data.instance.get_selected(true);
        let leaves = $.grep(objects, function (o) {return data.instance.is_leaf(o)});
        selectedItems = [];
        $.each(leaves, function (i, o) {
            selectedItems.push(o.text);
        });
        let selectedItemStr = selectedItems.join(",");

        // check if the previous filter is the same
        let tree_path = [userDir, "biopan", "selected_items.csv"].join('/');
        $.ajax({
            type: 'GET',
            url: 'get_csv_data.php', 
            data: {path: tree_path},
            dataType: 'text',
            success: function(result){
                res = result.split(",");
                res[res.length - 1] = res[res.length - 1].replace(/(\r\n|\n|\r)/gm,"");
                let equal = JSON.stringify(res) == JSON.stringify(selectedItems);
                if(!equal){
                    filterPathway(selectedItemStr);
                }else{
                    $('.graph-loading').hide();
                } 
            }
        });
    }


    //Update the filter on changed in the jstree
    let to = false;
    function changedFilter(data){
        $('.graph-loading').show();
        if(data.node && $("#search").val() != ""){ // if already a search (a filter in the search box)
            filterOnSearch(data);
        }else{
            if(to) {clearTimeout(to);}
            to = setTimeout(function(){
                defineFilter(data);
            }, 2000);
        }
    }


    //Update the filter on changed in the options
    function changedFilterInterface(){
        $(dtKey).hide();
        dtKey = '#' + pwLevel.value + '_' + subset.value;
        $(dtKey).jstree("check_all");
        $(dtKey).show();
        let value = $("#search").val();
        onSearch(value);
    }

    function showOptions(){
        $('#pw_level').show();
        $('#sub_set').show();
        $('#filter').show();
    }
    function hideOptions(){
        $('#pw_level').hide();
        $('#sub_set').hide();
        $('#filter').hide();
    }


    $("#class_reaction").on("changed.jstree", function (e, data) {
        if(dtKey == "#class_reaction"){
            if(to) {clearTimeout(to);}
            changedFilter(data);
        }
    });
    $("#species_reaction").on("changed.jstree", function (e, data) {
        if(dtKey == "#species_reaction"){
            changedFilter(data);
        }
    });
    $("#class_pathway").on("changed.jstree", function (e, data) {
        if(dtKey == "#class_pathway"){
            changedFilter(data);
        }
    });
    $("#species_pathway").on("changed.jstree", function (e, data) {
        if(dtKey == "#species_pathway"){
            changedFilter(data);
        }
    });
        
    
    //Function to relaunch the analysis if p-value and/or paired data options are not the default one (pValue=0.05, isPaired=No)
    function updateGraph(){
        $('.graph-loading').show();
        // call function to compute pathways and update graph
        $.ajax({
            type: 'GET',
            url: 'compute_pathway.php',
            data: {pw_level:pwLevel.value, p_value:pValue.value, is_paired:isPaired.value, userDir:userDir},
            dataType: 'text',
            success: function(result){
                $('.graph-loading').hide();
                openAnalysisType1(event, 'lipid-all', '<?php echo "$userDir"; ?>');
                if(result != null){
                    drawGraph();
                    onSearch("");
                    document.getElementById('search').value = "";
                }
            }
        });
    }


    //Highlight in green the active/suppressed nodes and edges
    function highLightSigPaths(highlighted_json_path){
        $.ajax({
            type: 'GET',
            url: 'get_json_data.php', 
            data: {path: highlighted_json_path},
            dataType: 'json',
            success: function(result){
                if(result){
                    let str_nodes = result['nodes'];
                    let str_edges = result['edges'];
                    let nodes = cy.$(str_nodes);
                    let edges = cy.$(str_edges);
                    edges.css({'overlay-color':'#68b0ab' ,'overlay-opacity':'0.4', 'overlay-padding':'10'});
                    nodes.css({'background-color':'#68b0ab', 'font-weight':'bold'});
                }
            }
        });
    }

    
    //Check that the 2 conditions to compare have the same number of samples
    function checkPairedGroups(diseaseGr, controlGr){
        let type_of_expr = document.getElementById('type_of_expr').value;
        if(type_of_expr){
            if(gr_freq[diseaseGr.value] != gr_freq[controlGr.value]){
                alert('Sample sizes in conditions must be equal!');
                return false;
            }
        }
        return true;
    }

    //Check that the 2 conditions to compare are different
    function checkSameGroups(diseaseGr, controlGr){
        if(diseaseGr.value == controlGr.value){
            alert('Conditions must be different!');
            return false;
        }
        return true;
    }


    //Export the graph in image or text format
    function exportGraph(fileType){
        let fileTitle = status(pwStatus.value)[1];
        let content;
        switch(fileType){
            case 'jpeg':
                content = cy.jpeg();
                break;
            case 'png':
                content = cy.png();
                break;
            case 'json':
            case 'txt':
                content = cy.json();
        }
        
        let fileName = 'ExportedGraph_' + fileTitle + '_' + diseaseGr.value + 'vs' + controlGr.value;
        let type;     
        let fileToSave;

        if(fileType == 'png'){
            saveAs(cy.png(), fileName);
        }else if(fileType == 'jpeg'){
            saveAs(cy.jpeg(), fileName);

        }else if(fileType == 'json'){
            // Create a blob of the data
                type = 'application/' + fileType;
                fileToSave = new Blob([JSON.stringify(content)], {
                    type: type,
                    name: fileName
                });
                saveAs(fileToSave, fileName);

        }else if(fileType == 'txt'){
            type = 'application/' + fileType;
            let output = 'Reactions chains' + '\t' + 'Z-score' + '\n';
            let source = '';
            let target = '';

            for(let i in content['elements']['edges']){
                for(let j in content['elements']['nodes']){
                    if(content['elements']['edges'][i]['data']['source'] == content['elements']['nodes'][j]['data']['id']){
                        source = content['elements']['nodes'][j]['data']['label'];
                    }
                    if(content['elements']['edges'][i]['data']['target'] == content['elements']['nodes'][j]['data']['id']){
                        target = content['elements']['nodes'][j]['data']['label'];
                    }
                }
                output += source + '->' + target + '\t' + content['elements']['edges'][i]['data']['weight'] + '\n';
            }
            fileToSave = new File([output], {
                type: 'text/plain',
                name: fileName
            });
            saveAs(fileToSave, fileName);
        }
    }


    let toJson = obj => obj.json();
    let toText = obj => obj.text();
    let tryPromise = fn => Promise.resolve().then(fn);
    
    let getStylesheet = name => {
        let convert = res => name.match(/[.]json$/) ? toJson(res) : toText(res);
        let data = fetch(`stylesheets/${name}`, {cache: 'no-cache'}).then( convert );
        return data;
    };
    
    let applyStylesheet = stylesheet => {
        if(typeof stylesheet === typeof ''){
            cy.style().fromString(stylesheet).update();
        }else {
            cy.style().fromJson(stylesheet).update();
        }
    };


    //On options changes
    pwLevel.addEventListener('change', function(){
        changedFilterInterface();
    });

    subset.addEventListener('change', function(){
        changedFilterInterface();
    });    
                
    pwType.addEventListener('change', function(){
        if(pwType.value == 'lp'){
            showOptions();
            changedFilterInterface();
        }else{
            hideOptions();
            drawGraph();
        }
    });

    pwStatus.addEventListener('change', function(){
        if(pwType.value == 'lp'){
            changedFilterInterface();
        }else{
            drawGraph();
        }
    });

    controlGr.addEventListener('change', function(){
        if(checkSameGroups(diseaseGr, controlGr) && checkPairedGroups(diseaseGr, controlGr)){
            changedFilterInterface();
        }
    });

    diseaseGr.addEventListener('change', function(){
        if(checkSameGroups(diseaseGr, controlGr) && checkPairedGroups(diseaseGr, controlGr)){
            changedFilterInterface();
        }
    });


    //On "Calculate pathways" button click
    let calcBtn = document.getElementById('calculate');
    calcBtn.addEventListener('click', function(){
        if(checkSameGroups(diseaseGr, controlGr) && checkPairedGroups(diseaseGr, controlGr)){
            updateGraph();
        }
    });
        
    
    //On "export" button click
    $('#export-btn').click(function() { 
        let file_type_to_save = getFileTypeToSave();
        let res = file_type_to_save[0]

        if(file_type_to_save[1] == 'graph'){ //export the graph
            exportGraph(res);
        }else{ //export the results tables
            let json_path = getSigPathwayTablePath(res);
            exportResultTables(json_path, res, subset.value, pwStatus.value, diseaseGr.value, controlGr.value, pValue.value);
        }
    });

  
    $('#modalEdgeInfor').on('hidden.bs.modal', function(){
        $(this).find('.modal-body select').empty();
    });


  //Change the displayed graph according to the clicked result table 
  // (e.g. display the "lipid subclass" graph after clicking on the "lipid sublclass results" table)
    $('#lp-cl').click(function() {
        if(pwType.value != 'lp' || pwLevel.value != 'class'){
            pwType.value = 'lp';
            pwLevel.value = 'class';
            showOptions();
            changedFilterInterface();
        }
    });
    $('#lp-sp').click(function() { 
        if(pwType.value != 'lp' || pwLevel.value != 'species'){
            pwType.value = 'lp';
            pwLevel.value = 'species';
            showOptions();
            changedFilterInterface();
        }
    });
    $('#all-fa').click(function() {
        if(pwType.value != 'fa'){
            pwType.value = 'fa';
            hideOptions();
            drawGraph();
        }
    });


    //Graph design
    let layouts = {
        coseBilkent:{
            name:'cose-bilkent',
            ready: undefined,
            stop: undefined,
            quality: 'default',
            fit: true,
            padding: 30,
            randomize: true,
            idealEdgeLength: 200,
        }
    };

    let prevLayout;
    let getLayout = name => Promise.resolve(layouts[name]);

    let applyLayout = layout => {
        if(prevLayout){
            prevLayout.stop();
        }
        let l = prevLayout = cy.makeLayout(layout);

        return l.run(function() {
            $("#loading").css('display','none');
        });
    }

    //Keep only the "cose-bilkent" layout. Keep functions to easily add another layout
    let applyLayoutFromSelect = () => Promise.resolve('coseBilkent').then(getLayout).then(applyLayout);
    //Keep only the "subclass.cycss" style. Keep functions to easily add another style
    let applyStyleSheetFromSelect = () => Promise.resolve('subclass.cycss').then(getStylesheet).then(applyStylesheet);

    // set style and layout
    tryPromise(applyStyleSheetFromSelect).then(applyLayoutFromSelect);

    let applyGraph = graphPath => {
        $.ajax({
            type: 'GET',
            url: 'get_json_data.php', 
            data: {path: graphPath},
            dataType: 'json',
            success: function(graph){
                if(graph != null){
                    cy.elements().remove();
                    cy.add(graph);
                    Promise.resolve('coseBilkent').then(getLayout).then(applyLayout);
                    
                    // highlight significant pathways if there are any
                    let highlighted_json_path = getHighLightedGraphPath();
                    highLightSigPaths(highlighted_json_path);

                    if(graph.nodes.length > 50){
                        // smaller text labels if more than 50 nodes
                        Promise.resolve('small.cycss').then(getStylesheet).then(applyStylesheet);
                    }else{
                        if(pwType.value == "fa" || pwLevel.value == "species"){
                            Promise.resolve('species.cycss').then(getStylesheet).then(applyStylesheet);
                        }else{
                            Promise.resolve('subclass.cycss').then(getStylesheet).then(applyStylesheet);
                        }
                    }
                }else{
                    cy.elements().remove();
                }
            },
            error: function(){
                cy.elements().remove();
            }
        });
    }

    
    //Draw the requested graph
    function drawGraph(){
        Promise.resolve(userDir).then(getGraphDatasetPath).then(applyGraph);
        let highlighted_json_path = getHighLightedGraphPath();
        highLightSigPaths(highlighted_json_path);
        updateTables(subset.value, pwStatus.value, diseaseGr.value, controlGr.value);
    }

    
    //Create HTML to display the selected/non-selected lipid molecular species after clicking on a edge (reaction)
    function addToListBox(elem, nonselected, selected){
        if(elem != null){
            if(nonselected.length > 0){  
                $.each(nonselected, function(i, e) {
                    elem.append($('<option></option>')
                    .attr({'value': e}).text(e));
                });             
            }
            if(selected.length > 0){
                $.each(selected, function(i, e) {
                    elem.append($('<option></option>')
                    .attr({'value': e, 'selected': true}).text(e)); 
                });
            }
        }
    }

    
    if(cy == null){
        // init graph
        initGraph();
    }

    //Initialisation of the graph
    function initGraph(){
        cy = cytoscape({
            container: elem,
            elements: [],
            ready: function (e) {
                let cy = e.cy;
                // show a link when user clicks on the node
                cy.on('tap', 'node', function(){
                    let label = this.data('label');
                    let species = '';

                    let url_1 = 'https://lipidmaps.org/data/structure/LMSDFuzzySearch.php?Name=';
                    let url_2 = '&SortResultsBy=Name';
                    
                    if(this.data('lm_id').length !== 0){ // for nodes of subclasses graphs that have a LMID
                        url_1 = 'https://lipidmaps.org/data/structure/LMSDSearch.php?Mode=ProcessClassSearch&LMID=';
                        url_2 = '';
                        species = this.data('lm_id');

                    }else if(label.match('O-|P-') !== null){ // reconstruct nomenclature for plasmalogens (O-PC(34:0) -> PC(O-34:0))
                        let plasmalogens = 'O-';
                        if(label.match('P-')){
                            plasmalogens = 'P-';
                        }
                        label = label.replace(plasmalogens, '');
                        species = label.split('(');
                        if(pwType.value == 'lp' && pwLevel.value == 'class'){ // subclasses graphs: there is no structure
                            species = species[0].concat('(', plasmalogens);
                        }else{
                            species = species[0].concat('(', plasmalogens, species[1]);
                        }
                        
                    
                    }else if(label.match('Cer|dhCer|SM|dhSM') !== null){ // for cer/dhcer/sm/dhsm
                        let backbone = 'd18:1';
                        if(label.match('dh')){
                            backbone = 'd18:0';
                        }
                        let lp = 'Cer';
                        if(label.match('SM')){
                            lp = 'SM';
                        }
                        species = label.split('(');
                        species = lp.concat('%28', backbone, '%2F', species[1]);
                        species.replace(':', '%3A').replace(')', '%29');
                        url_1 = 'https://lipidmaps.org/search/quicksearch.php?Name='; // use another link to get more specifics results
                        url_2 = '';
                    }else{
                        species = label;
                    }
                    let url = url_1.concat(species, url_2);
                    
                    try { // your browser may block popups
                        window.open(url);
                    } catch(e){ // fall back on url change
                        window.location.href = url;
                    }
                });                    
            }
        });
                      
        cy.panzoom();

        $('#cy').css('overflow', 'hidden');
        cy.on('mouseover', 'node', function(e){
            let sel = e.target;
            cy.elements().difference(sel.outgoers()).not(sel).addClass('semitransp');
            sel.addClass('highlight').outgoers().addClass('highlight');
            $('#cy').css('cursor', 'pointer');

            if(pwType.value == 'lp' && pwLevel.value == 'class'){
                let makeTippy = function(node, text) {
                return tippy(node.popperRef(), {
                    html: (function() {
                        let div = document.createElement('div');
                        div.innerHTML = text;
                        return div;
                    })(),
                    trigger: 'manual',
                    arrow: true,
                    placement: 'bottom',
                    hideOnClick: false,
                    multiple: true,
                    sticky: true
                    }).tooltips[0];
                };
                tippyNode = makeTippy(this, this.data('name'));    
                tippyNode.show();
            }
        });
                    
        cy.on('mouseout', 'node', function(e){
            let sel = e.target;
            cy.elements().removeClass('semitransp');
            sel.removeClass('highlight').outgoers().removeClass('highlight');
            $('#cy').css('cursor', 'default');
            if(tippyNode != null){
                tippyNode.hide();
            }
        });
                
        cy.on('mouseover', 'edge', function(e){
            $('#cy').css('cursor', 'pointer');
            e.target.css({'width': '6px'});
        });
                    
        cy.on('mouseout', 'edge', function(e){
            $('#cy').css('cursor', 'default');
            e.target.css({'width': '4px'});
            
        });
                    
        // show selected/non-selected lipid molecular species after clicking on an edge
        cy.on('tap', 'edge', function(){                 
            re = this.data('source');
            pro = this.data('target');
            let tokens = userDir.split('/');
            let userId = tokens.pop();
            let react = [re,pro].join(',').toLowerCase();
            let json_path = [userDir, 'biopan', react.concat('.json')].join('/');
                                    
            return $.ajax({
                type: 'GET',
                url: 'get_json_data.php', 
                data: { path: json_path },
                dataType: 'json',
                success: function(result){             
                    if(result != null){
                        listBoxRe = $('#dual-list-box-re').bootstrapDualListbox({
                            nonSelectedListLabel: 'Non-selected',
                            selectedListLabel: 'Selected',
                            preserveSelectionOnMove: 'moved',
                            moveOnSelect: false,
                            nonSelectedFilter: ''
                        });

                        listBoxPro = $('#dual-list-box-pro').bootstrapDualListbox({
                            nonSelectedListLabel: 'Non-selected',
                            selectedListLabel: 'Selected',
                            preserveSelectionOnMove: 'moved',
                            moveOnSelect: false,
                            nonSelectedFilter: ''
                        });

                        let non_selected_re = result['non_selected_re'];
                        
                        let non_selected_re_arr = [];
                        if(non_selected_re.length > 0){
                            non_selected_re_arr = non_selected_re.split(',');
                        }
                                                
                        let selected_re = result['selected_re'];
                        let selected_re_arr = [];
                        if(selected_re.length > 0){
                            selected_re_arr = selected_re.split(',');
                        }    
                        all_re_arr = non_selected_re_arr.concat(selected_re_arr);
                        let non_selected_pro= result['non_selected_pro'];   
                        let non_selected_pro_arr = [];
                        if(non_selected_pro.length > 0){
                            non_selected_pro_arr = non_selected_pro.split(',');
                        }
                                                
                        let selected_pro = result['selected_pro'];
                        let selected_pro_arr = [];
                        if(selected_pro.length > 0){
                            selected_pro_arr = selected_pro.split(',');
                        }   
                        all_pro_arr = non_selected_pro_arr.concat(selected_pro_arr);    
                                                
    
                        addToListBox(listBoxRe, non_selected_re_arr, selected_re_arr);
                        
                        addToListBox(listBoxPro, non_selected_pro_arr, selected_pro_arr); 
                                                
                        listBoxRe.bootstrapDualListbox('refresh');
                        listBoxPro.bootstrapDualListbox('refresh');
                                            
                        $('#modalEdgeInfor').modal();
                        $('.filter.form-control').hide();
                        $('.btn-group.buttons').hide();  
                    }
                },
                error: function(xhr, status, err) { alert('Ajax request failed'); }
            });
        });     
    }
        
    // resizing    
    cy.on('ready', function () {
        updateBounds();
    });

    let updateBounds = function () {
        let bounds = cy.elements().boundingBox();
        $('#cyContainer').css('height', bounds.h + 300);
        cy.center();
        cy.resize();
        // fix the Edgehandles
        $('#cy').cytoscapeEdgehandles('resize');
    };


    drawGraph();
    $(dtKey).show();
    openAnalysisType1(event, 'lipid-all', '<?php echo "$userDir"; ?>');
}
