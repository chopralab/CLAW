/* To check the format fie */
function fileValidation() { 
    var fileInput = document.getElementById('dataFile');
    var filePath = fileInput.value;
    var allowedExtensions = /(\.csv)$/i;
    if (!allowedExtensions.exec(filePath)) {
        alert('Invalid file type: BioPAN requires a CSV file');
        fileInput.value = '';
        return false; 
    }  
} 


/* To download the demo file */
function download(url, filename) {
    fetch(url).then(function(t) {
        return t.blob().then((b)=>{
            var a = document.createElement("a");
            a.href = URL.createObjectURL(b);
            a.setAttribute("download", filename);
            a.click();
        }
        );
    });
}