from jinja2 import Environment
import os

base = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>NblockMatcher</title>
  <script src="https://cdn.jsdelivr.net/npm/igv@2.15.10/dist/igv.min.js"></script>
</head>
<body>
<div>
    <h1>NblockMatcher result visulaziation</h1>
    <div>
        <div id="igv-div"></div>
    </div>
</div>
</body> 
<footer>
<script>
var igvDiv = document.getElementById("igv-div");

var options =
{
    "reference": {
        "name": "genome",
        "fastaURL": "{{genome_fasta}}",
        "indexed": false,
        "showIdeogram": false,
    },
    "showChromosomeWidget": false,
    "loadDefaultGenomes": false,
    "tracks": [
    {
        "name": "NB_sites",
        "type": "annotation",
        "format": "gff3",
        "displayMode": "expanded",
        "height": 120,
        "url": "{{prefix}}.gff",
        "order": 2
    },
    {
        "name": "NB_arcs",
        "url": "{{prefix}}.interact",
        "type": "interaction",
        "format": "interact",
        "arcType": "proportional",
        "height": 80,
        "showBlocks": true,
        "order": 1
    },
    ]
};

igv.createBrowser(igvDiv, options)
        .then(function (browser) {
             igv.browser = browser;
        })

</script>

</footer>
</html> 
"""


def to_html(prefix, fasta):

    template = Environment().from_string(base)

    with open(prefix + '.html', 'w') as f:
        f.write(template.render(
            prefix=os.path.basename(prefix),
            genome_fasta=os.path.relpath(fasta, os.path.dirname(prefix)),
        ))
