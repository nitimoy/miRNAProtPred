// Function to validate the sequence field
function validateSequence() {
    var sequenceInput = document.getElementById('sequence');
    var sequenceValue = sequenceInput.value;
    
    // Use a regular expression to check if the input contains any non-alphabet characters except spaces
    var containsNonAlphabet = /[^a-zA-Z\s]/.test(sequenceValue);
    
    if (containsNonAlphabet) {
        // If it contains non-alphabet characters (other than spaces), show an error message
        alert('Error: Sequence should only contain alphabets and spaces.');
        sequenceInput.focus(); // Put focus back on the sequence field
        return false; // Prevent form submission
    }
    
    return true; // Allow form submission if only alphabets and spaces are found
}


// Function to show the waiting modal
function showWaitingModal() {
    if (validateSequence()) {
        $('#waitingModal').modal({
            backdrop: 'static',
            keyboard: false
        });
        return true; // Ensure form submission
    }
    return false; // Prevent form submission if validation fails
}
