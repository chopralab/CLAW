// To launch BioPAN
function submitForm() {  
    var samples = [];
    var groups = []; 
    
    $("input[name*='samples']").each(function() { samples.push( $(this).val());});
    $("input[name*='groups']").each(function() { groups.push( $(this).val());});
    
    var str_groups = groups.join(",");
    var str_samples = samples.join(",");
    var userID = document.getElementById("userid").value;
    
    if(groups.length < 2){
        $(".less-than-one-sample-alert").show();
        return;
    }

    // To display progress bar
    $(".overlay").show();
    $(".loading").css('top','54%');
    $(".loading").css('left','47%');
    
    // Transfering form information to different page without page refresh
    $.ajax({
        type: "POST",
        url: "processing.php", 
        data: {userID: userID, samples: str_samples, groups: str_groups},
        dataType: "text",
        success: function(res){
            $(".overlay").hide();
            window.location = res;
        },
        error:function(jqXHR, textStatus, errorThrown){
            alert("Error type" + textStatus + "occured, with value " + errorThrown);
        }
    });
}

                
//Assign group: to define which group contain at least 2 samples
function getColor(){
    let values = {};
    for(i = 0; i < valuesTable.length; i++){
        if(!(valuesTable[i].value in values)){
            values[valuesTable[i].value] = 1;
        }else{
            values[valuesTable[i].value] = values[valuesTable[i].value] + 1;
        }
    }

    let valid = [];
    for(i = 0; i < valuesTable.length; i++){
        let id = 'input#valuesTable'.concat(i);
        if(values[valuesTable[i].value] < 2){
            $(id).css('background-color','#feb9b9');
        }else{
            $(id).css('background-color','#c7feb9');
            if(valid.indexOf(valuesTable[i].value) === -1){
                valid.push(valuesTable[i].value);
            }
        }
    }
    if(valid.length < 2){
        document.getElementById('assignGroups').disabled = true;
    }else{
        document.getElementById('assignGroups').disabled = false;
    }
}