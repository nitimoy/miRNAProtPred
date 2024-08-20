$(document).ready(function() {
    var currentPath = window.location.pathname;
    $('#navbarNav .nav-link').each(function() {
        var linkPath = $(this).attr('href');
        if (linkPath === currentPath || linkPath === currentPath + '/') {
            $(this).addClass('active');
        }
    });
});
