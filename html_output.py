from jinja2 import Environment
import os
import re
import shutil
# todo update igv when fixed version is released
base = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>NblockMatcher</title>
  <script src="https://cdn.jsdelivr.net/npm/igv@2.15.10/dist/igv.min.js"></script>
  <style>
    .clickme {color: blue; text-decoration: underline;}
    .clickme:hover { cursor: pointer; }

    .scroll {
      overflow-y: auto; 
      height: 500px; 
    }
    
    .scroll thead th {
      position: sticky; 
      top: 0px; 
    }
    
    th {
      background: #f3f3f3;
    }
    
  </style>
</head>
<body>
<div>
    <h1>NblockMatcher result visualisation</h1>
    <div>
        <div id="igv-div"></div>
    </div>
    <div class="scroll">
        {{dfhtml}}
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
    {% if gff != '' %}
        {
            "name": "gff",
            "type": "annotation",
            "format": "gff3",
            "displayMode": "expanded",
            "url": "{{prefix}}_spec.gff",
            "order": 3
        },
    {% endif %}
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


def to_html(prefix, fasta, gff='', df=None):
    if gff != '':
        shutil.copyfile(gff, f'{prefix}_spec.gff')

    def func2(x, offset=200):
        sinds = [int(i) for i in re.findall(r'\d+', x.sites)]
        mi = max(1, min(sinds) - offset)
        ma = max(sinds) + offset

        return (
            f"""<tr><td><span class="clickme" onClick="igv.browser.search('{x.sequence_name}:{mi}-{ma}')">{round(x.score, 2):0.2f}</span></td>"""
            f"<td>{x.sequence_name}</td>"
            f"<td>{x.orientation}</td>"
            f"<td>{x.matched_sequences}</td>"
            f"<td>{min(sinds)}-{max(sinds)}<td></tr>"
        )

    if df is None:
        dfhtml = ''
    else:
        dfhtml = ('<table>'
                  '<thead><tr><th>score</th><th>seq</th><th>orientation</th><th>match</th><th>range</th></tr></thead>'
                  '<tbody>') + '\n'.join([func2(i) for i in df.itertuples()]) + '\n</tbody><table>'

    template = Environment().from_string(base)

    with open(prefix + '.html', 'w') as f:
        f.write(template.render(
            prefix=os.path.basename(prefix),
            genome_fasta=os.path.relpath(fasta, os.path.dirname(prefix)),
            gff=gff,
            dfhtml=dfhtml,
        ))
