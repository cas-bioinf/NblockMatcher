import re


def site2minmax(site):
    x = [int(i) for i in re.findall(r'\d+', site)]
    return min(x), max(x)