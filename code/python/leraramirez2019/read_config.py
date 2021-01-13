#!/usr/bin/env python2
# Changed by manu to work on python2!

"""
    Read cytosim configuration files

Syntax:
   
    read_config.py FILE
    read_config.py - FILE
    
Description:

    Reads a config file and prints a formatted copy to standard output.
    Two formats are possible: vertical (default) or horizontal (second syntax)

F. Nedelec 11.2013--2015
"""

import sys, os, io

err = sys.stderr
out = sys.stdout


def format_value(val):

    if isinstance(val, str):
        return val
    if isinstance(val, dict):
        r = ""
        s = "{ "
        for k in sorted(val):
            r += s + str(k) + " = " + format_value(val[k])
            s = "; "
        return r + "; }"
    try:
        r = ""
        s = ""
        for k in val:
            r += s + format_value(k)
            s = ", "
        return r
    except:
        return repr(val)


class Instruction:
    keys = []
    cnt  = 1
    pam  = {}
    def __init__(self, key):
        self.keys = []
        self.keys.append(key)
        self.keys.append('')
        self.keys.append('')
        self.field = []
    def __repr__(self):
        return self.horizontal()
    def horizontal(self):
        res = self.keys[0] + " "
        if self.cnt != 1:
            res += repr(self.cnt) + " "
        res += self.keys[1] + " " + self.keys[2] + " { "
        for p in sorted(self.pam):
            res += str(p) + " = " + format_value(self.pam[p]) + "; "
        res += "}"
        return res
    def vertical(self):
        res = self.keys[0] + " "
        if self.cnt != 1:
            res += repr(self.cnt) + " "
        res += self.keys[1] + " " + self.keys[2] + "\n{\n"
        for p in sorted(self.pam):
            res += "    " + str(p) + " = " + format_value(self.pam[p]) + ";\n"
        res += "}\n"
        return res
    def value(self, key, index=0):
        if key in self.pam:
            val = self.pam[key]
            if isinstance(val, list):
                return val[index]
            else:
                return val
        return ""
    def values(self):
        return self.pam;


def get_character(fid):
    # THE ONLY CHANGE FOR PYTHON2 to work
    return str(fid.read(1))


def get_hexadecimal(fid, s):
    res = s
    pos = fid.tell()
    c = get_character(fid)
    while c and c in "0123456789ABCDEFabcdef":
        res += c
        pos = fid.tell()
        c = get_character(fid)
    fid.seek(pos, 0)
    return res


def get_number(fid, s):
    """
    Read a number with decimal point and optional exponent
    """
    res = s
    pos = fid.tell()
    c = get_character(fid)
    if c == 'x':
        return get_hexadecimal(fid, s+c)
    if c == '-' or c == '+':
        res += c
        pos = fid.tell()
        c = get_character(fid)
    no_point = 1
    while c.isdigit() or ( c=='.' and no_point ):
        if c=='.':
            no_point = 0
        res += c
        pos = fid.tell()
        c = get_character(fid)
    if c == 'e' or c == 'E':
        return res + get_number(fid, 'e')
    if c == 'P':
        return res
    fid.seek(pos, 0)
    return res


def delimiter(c):
    if c == '(': return ')'
    if c == '{': return '}'
    if c == '[': return ']'
    if c == '"': return '"'
    return 0


def get_until(fid, e):
    res = ''
    c = get_character(fid)
    while c and c != e:
        res += c
        c = get_character(fid)
    return res


def get_block(fid, s, e):
    """
    Return a block including the enclosing delimiters
    """
    res = s
    c = get_character(fid)
    while c:
        if c == e:
            res += c
            return res
        if delimiter(c):
            res += get_block(fid, c, delimiter(c))
        else:
            res += c
        c = get_character(fid)
    err.write("  Error: unclosed block\n")
    sys.exit(1)
    return res


def valid_token_char(c):
    return c.isalnum() or c == '_' or c == '.'



def get_token(fid):
    """
    Extract the next token from the file
    """
    c = get_character(fid)
    while c.isspace() and c != '\n':
        c = get_character(fid)
    if c == '%':
        fid.readline()
        return "%"
    if delimiter(c):
        return get_block(fid, c, delimiter(c))
    if c.isdigit() or c == '-' or c == '+':
        return get_number(fid, c)
    res = c
    if valid_token_char(c):
        while c:
            pos = fid.tell()
            c = get_character(fid)
            if not valid_token_char(c):
                fid.seek(pos, 0)
                break
            res += c
    #print("token "+ res)
    return res


def uncode(arg):
    if not isinstance(arg, str):
        try:
            return arg.decode('utf-8')
        except:
            pass
    return arg


def simplify(arg):
    """
    Simplify the values of imbricated sets and dictionaries
    """

    val = uncode(arg)
    if isinstance(val, str):
        try:
            return int(val)
        except:
            pass
        try:
            return float(val)
        except:
            pass
        s = val.strip()
        #remove outer delimiters:
        if len(s) > 3 and delimiter(s[0]) == s[-1]:
            s = s[1:-1].strip()
        return s
    if isinstance(val, dict):
        r = {}
        for k in val:
            r[k] = simplify(val[k])
        return r
    try:
        if len(val) == 1:
            return simplify(val[0])
        r = []
        for a in val:
            r.append(simplify(a))
        return r
    except:
        return val


def file_object(arg):

    #print("new_file_object", type(arg))
    try:
        return io.StringIO(arg)
    except:
        return io.StringIO(unicode(arg))


def read_list(fid):
    """
    Process a list of key=values
    """
    dic = {}
    key = ''
    val = ''
    tak = ''
    tok = get_token(fid)
    while tok:
        #print('key ', tok)
        if tok == '=':
            key = tak;
            dic[key] = []
            val = ''
        elif tok == ';' or tok == '\n' or tok[0] == '%':
            if key:
                dic[key].append(val)
            key = ''
            val = ''
        elif tok == ',':
            if key:
                dic[key].append(val)
            val = ''
        elif delimiter(tok[0]):

            val = read_list(file_object(tok[1:-1]))
            if not val:
                val = tok
        else:
            if key:
                if val:
                    val += ' ' + tok
                else:
                    val = tok
            else:
                tak = tok
        tok = get_token(fid)
    if key and val:
        dic[key].append(val)
    return dic



def parse_config(fid):
    """
    return `pile` resulting from parsing the file
    """
    lev = 0
    pile = []
    while fid:
        tok = get_token(fid)
        if not tok:
            break
        #print('level', lev, '  token ', tok)
        if tok[0] == '%':
            pass
        elif tok == '\n':
            pass
        elif tok in [ 'set', 'new', 'run', 'call', 'delete', 'change', 'cut', 'mark' ]:
            if lev == 3 and 'cur' in locals():
                cur.pam = {}
                pile.append(cur)
            cur = Instruction(tok)
            lev = 1
        elif tok == 'repeat':
            cnt = get_token(fid)
            blok = get_token(fid)
            sfid = io.StringIO(blok[1:-1])
            pile.extend(parse_config(sfid))
        elif lev == 1 and tok.isdigit():
            cur.cnt = int(tok)
        elif lev == 1 and tok[0] == '[':
            cur.cnt = tok
        elif lev == 1:
            cur.keys[lev] = tok
            lev = 2
        elif lev == 2 and tok.isdigit():
            cur.inx = int(tok)
        elif lev == 2:
            if tok == ':':
                cur.field = get_token(fid);
            else:
                cur.keys[lev] = tok
                lev = 3
        elif lev == 3 and delimiter(tok[0]):
            if cur.field:
                cur.pam = {cur.field:tok}
            else:
                pam = read_list(file_object(tok[1:-1]))
                cur.pam = simplify(pam)
                #print("simplified\n", pam, "\n", cur.pam)
            pile.append(cur)
            lev = 0
    return pile



def parse(arg):
    """
    Returns a dictionnary obtained by parsing the file at `arg`
    """
    if isinstance(arg, str):
        f = open(arg, 'r')
        pile = parse_config(f)
        f.close()
    else:
        pile = parse_config(arg)
    return pile

#------------------------------------------------------------------------

def format_horizontal(pile, prefix=''):
    """
    Returns a compact representation of the hierachical dictionnary
    """
    res = ''
    for i in pile:
        res += prefix + i.horizontal() + "\n"
    return res


def format_vertical(pile, prefix=''):
    """
    Returns a neat representation of the hierachical dictionnary
    """
    res = ''
    for i in pile:
        res += prefix + i.vertical() + "\n"
    return res

#------------------------------------------------------------------------

def get_command(pile, keys):
    """
        Return all commands corresponding to = [ command, class, name ]
        The keys can be '*' to specify a wildcard matching
        """
    if len(keys) != 3:
        return "Error: get_command(arg, keys) expects `keys` to be of size 3"
    res = []
    for p in pile:
        if keys[0]=='*' or p.keys[0]==keys[0]:
            if keys[1]=='*' or p.keys[1]==keys[1]:
                if keys[2]=='*' or p.keys[2]==keys[2]:
                    res.append(p)
    if len(res) == 1:
        return res[0]
    else:
        return res


def get_value(arg, keys):
    """
    Return the value speficied by keys = [ command, class, name, parameter ]
    The first 3 keys can be '*' for wildcard matching
    """
    if isinstance(arg, str):
        pile = parse(arg)
    else:
        pile = arg
    if len(keys) != 4:
        return "Error: get_value(arg, keys) expects `keys` to be of size 4"
    for p in pile:
        if keys[0]=='*' or p.keys[0]==keys[0]:
            if keys[1]=='*' or p.keys[1]==keys[1]:
                if keys[2]=='*' or p.keys[2]==keys[2]:
                    try:
                        return p.pam[keys[3]]
                    except KeyError:
                        if keys[3]=='*':
                            return p.pam
    return "Unspecified"

def get_vals(pile,extract):
    out = list()
    for i in range(0, len(extract), 3):
        com = extract[i]
        what = extract[i + 1]
        if len(com)==3:
            obj = get_command(pile, com)
        else:
            obj = get_command(pile, com[:3])
            obj = obj[com[-1]]
        if len(what)==1:
            out.append(obj.value(*what))
        elif len(what)>1:
            out.append(obj[what[1]].value(what[0]))
        else:
            out.append(obj.cnt)
    return out

#------------------------------------------------------------------------

def main(args):
    files = []
    horizontal = False
    
    for arg in args:
        if os.path.isfile(arg):
            files.append(arg)
        elif arg=="-":
            horizontal = True
        else:
            err.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    
    if not files:
        err.write("  Error: you must specify a file to read\n")
        sys.exit()

    for f in files:
        out.write("\n")
        pile = parse(f)
        if horizontal:
             out.write(format_horizontal(pile))
        else:
             out.write(format_vertical(pile))
        #print(f.rjust(40) + " " + value(pile, ['set', 'space', '*', 'geometry']) )


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])
