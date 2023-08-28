import signal

from jinja2 import Environment
import os
import re
import shutil
import gzip
from src.utils import site2minmax

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
      font-family: monospace;
    }
    
    .scroll thead th {
      position: sticky; 
      top: 0px; 
    }
    
    th {
      background: #f3f3f3;
    }
    
    .highlight {
     background: #9e9e9e;
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
var rows = document.querySelectorAll('#matches tr');

function highlightMe(ele) {
      rows.forEach(function(r) {
        r.classList.remove('highlight');
      })
      ele.parentNode.parentNode.classList.add('highlight');
    }
</script>

<script>
var igvDiv = document.getElementById("igv-div");

var options =
{
    "reference": {
        "name": "genome",
        "fastaURL": "{{prefix}}_genome.fa",
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
            "order": 3,
            "filterTypes": {{filterTypes}},
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


def to_html(prefix, fasta, gff='', df=None, filter_types=('region',)):
    if fasta.endswith('.gz'):
        with gzip.open(fasta, 'rt') as f, open(f'{prefix}_genome.fa', 'w') as fo:
            fo.write(f.read())
    else:
        shutil.copyfile(fasta, f'{prefix}_genome.fa')

    if gff != '':
        if gff.endswith('.gz'):
            with gzip.open(gff, 'rt') as f, open(f'{prefix}_spec.gff', 'w') as fo:
                fo.write(f.read())
        else:
            shutil.copyfile(gff, f'{prefix}_spec.gff')

    def func2(x, offset=200):
        mi, ma = site2minmax(x.sites)
        mi_s = max(1, mi - offset)
        ma_s = ma + offset

        return (
            f"""<tr><td><span class="clickme" onClick="highlightMe(this); igv.browser.search('{x.sequence_name}:{mi_s}-{ma_s}')">{round(x.score, 2):0.2f}</span></td>"""
            f"<td>{x.sequence_name}</td>"
            f"<td>{x.orientation}</td>"
            f"<td>{x.matched_sequences}</td>"
            f"<td>{mi}-{ma}</td>"
            f"<td>{x.dists}</td></tr>"
        )

    if df is None:
        dfhtml = ''
    else:
        dfhtml = ('<table id="matches">'
                  '<thead><tr><th>score</th><th>seq</th><th>orientation</th><th>match</th><th>range</th><th>dists</th></tr></thead>'
                  '<tbody>') + '\n'.join([func2(i) for i in df.itertuples()]) + '\n</tbody><table>'

    template = Environment().from_string(base)

    with open(prefix + '.html', 'w') as f:
        f.write(template.render(
            prefix=os.path.basename(prefix),
            gff=gff,
            dfhtml=dfhtml,
            filterTypes=str(list(filter_types)),
        ))
