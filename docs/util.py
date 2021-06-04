def pp(x):
    '''
    Hacky way to pretty-print nested dict __repr__'s
    '''
    if not isinstance(x, str):
        x = x.__repr__()
    delims = {"'" : "'", '"' : '"', "(" : ")", "[" : "]"}
    out = []
    indent = [0]
    inside = []
    prefix = 0
    name = True
    for a in x:
        if prefix == 0 and a == " ":
            continue
        if a == "}":
            prefix = 0
            out.append("\n" + " " * (sum(indent) - 1))
            indent = indent[:-1]
        else:
            prefix += 1
        out.append(a)
        if a == '{':
            indent.append(min(prefix, 4))
            prefix = 0
            out.append("\n" + " " * sum(indent))
            inside = []
        if len(inside) > 0 and inside[-1] in delims:
            if a == delims[inside[-1]]:
                inside = inside[:-1]
        else:
            if a in delims:
                inside.append(a)
            if name:
                if a == ':':
                    name = False
            else:
                if a == ',':
                    name = True
                    prefix = 0
                    out.append("\n" + " " * sum(indent))
    out = "".join(out)
    print(out)
