<html>
<head>
<title>Stdlib and Utilities</title>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='StdlibAndUtilities.php'>

<h2>Stdlib and Utilities</h2>

<code>PythiaStdlib</code> and <code>PythiaComplex</code> only exist as 
header files, collecting some simple declarations and inline utilities.

<p/>
<code>PythiaStdlib</code> collects the <code>include</code> and 
<code>using</code> statements that are required by most other classes 
to access the C++ <code>Stdlib</code> containers and methods, such as
<code>string</code>, <code>vector</code>, <code>map</code>, some 
mathematical functions, and input/output streams and formats.
It defines <code>M_PI</code> if this is not already done. 

<p/>
There are also a few inlined functions: <code>pow2(x)</code>, 
<code>pow3(x)</code>, <code>pow4(x)</code> and <code>pow5(x)</code> for 
small integer powers, and <code>sqrtpos(x)</code> where a 
<code>max(0., x)</code> ensures that one does not take the square root 
of a negative number.

<p/>
<code>PythiaComplex</code> defines a <code>complex</code> data type 
by a <code>typedef std::complex&lt;double&gt;</code>.

</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
