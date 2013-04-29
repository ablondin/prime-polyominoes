#!/usr/bin/python
import sys
from prime import *

# For testing
DRAGON = [(1,0), (2,0), (3,0), (5,0), (6,0), (3,1), (6,1), (0,2), (1,2), (3,2),\
          (4,2), (4,4), (5,2), (6,2), (0,3), (6,3), (7,3), (0,4), (1,4), (2,4),\
          (3,4), (5,4), (6,4), (7,4), (3,5), (6,5), (1,6), (3,6), (5,6), (6,6),\
          (1,7), (2,7), (3,7)]
 
# -------------------- #
# Handling polyominoes #
# -------------------- #

def canonical(p):
    mx = min(map(lambda v: v[0], p))
    my = min(map(lambda v: v[1], p))
    return sorted(map(lambda v: (v[0]-mx, v[1]-my), p))

def rotate(p):
    return canonical(map(lambda v: (v[1], -v[0]), p))

def reflect(p):
    return canonical(map(lambda v: (v[0], -v[1]), p))

def expand(p):
    result = []
    for (x,y) in p:
        for (dx,dy) in ((-1,0),(1,0),(0,-1),(0,1)):
            if p.count((x+dx,y+dy)) == 0:
              result.append(canonical(p + [(x+dx,y+dy)]))
    return result

def neighbors(cell):
    (x,y) = cell
    return [(x+1,y),(x,y+1),(x-1,y),(x,y-1)]

def has_hole(polyomino):
    w = max(map(lambda c: c[0], polyomino))
    h = max(map(lambda c: c[1], polyomino))
    complement = [(x,y) for x in xrange(-1, w+2) for y in xrange(-1,h+2) if (x,y) not in polyomino]
    queue = [complement[0]]
    visited = set([])
    while queue:
        cell = queue.pop()
        visited.add(cell)
        for neighbor in neighbors(cell):
            if neighbor in complement and neighbor not in visited:
                queue.append(neighbor)
    return len(visited) != len(complement)

def polyomino_to_boundary_word(polyomino):
    w = max(map(lambda c: c[0], polyomino))
    h = max(map(lambda c: c[1], polyomino))
    (x,y) = min(polyomino, key=lambda c: (c[1] != 0, c[0]))
    # North-east
    q1, x, y = quarter_word(polyomino, w, h, x, y, 0)
    q2, x, y = quarter_word(polyomino, w, h, x, y, 1)
    q3, x, y = quarter_word(polyomino, w, h, x, y, 2)
    q4, x, y = quarter_word(polyomino, w, h, x, y, 3)
    return q1 + q2 + q3 + q4

def neighbour((x,y), s):
    if   s == 0: return (x + 1, y)
    elif s == 1: return (x, y + 1)
    elif s == 2: return (x - 1, y)
    elif s == 3: return (x, y - 1)

def quarter_word(polyomino, w, h, x, y, s):
    word = [s]
    if   s == 0: condition = lambda (x,y,w,h): x != w
    elif s == 1: condition = lambda (x,y,w,h): y != h
    elif s == 2: condition = lambda (x,y,w,h): x != 0
    elif s == 3: condition = lambda (x,y,w,h): y != 0
    current = s
    while condition((x,y,w,h)) or current != s:
        left = (current + 1) % 4
        right = (current - 1) % 4
        if neighbour((x,y), current) in polyomino:
            (x,y) = neighbour((x,y), current)
            if neighbour((x,y), right) in polyomino:
                (x,y) = neighbour((x,y), right)
                current = right
        else:
            current = left
        word.append(current)
    return word, x, y

# --------- #
# Iterators #
# --------- #

def polyomino_iterator(max_cells):
    yield ('0123',[(0,0)])
    polyominoes = [[(0,0)]]
    for _ in xrange(0, max_cells - 1):
        new_polyominoes = []
        for p in polyominoes:
            for q in expand(p):
                dup = False
                for r in xrange(8):
                    if r == 4:
                        q = reflect(q)
                    if q in new_polyominoes:
                        dup = True
                        break
                    q = rotate(q)
                if not dup and not has_hole(q):
                    yield (''.join(map(str, polyomino_to_boundary_word(q))),q)
                    new_polyominoes.append(q)
        polyominoes = new_polyominoes

def polyominoes_iterator_from_file(input_filename):
    import ast
    input = open(input_filename, 'r')
    for line in input:
        (word,polyomino) = line.strip().split(':')
        word = tuple(map(int, word))
        polyomino = ast.literal_eval(polyomino)
        if is_simple(word[:-1]):
            bw = BoundaryWord(word)
            yield (word,polyomino)
    input.close()

def composed_polyominoes_iterator_from_file(input_filename):
    for (word,polyomino) in polyominoes_iterator_from_file(input_filename):
        bw = BoundaryWord(word) 
        if not bw.is_prime():
            yield (word,polyomino)

# -------------- #
# Output to file #
# -------------- #

def write_free_polyominoes_to_file(max_cells, output_filename=None, sorted=True):
    if output_filename is None:
        output_filename = 'free-polyomino-%s.txt'%max_cells
    f = open(output_filename, 'w')
    polyominoes = list(polyomino_iterator(max_cells))
    if sorted:
        polyominoes.sort(key=lambda p: (len(p[1]),len(p[0])))
    for (boundary_word,polyomino) in polyominoes:
        f.write('%s:%s\n'%(boundary_word,polyomino))
    f.close()

def write_composed_polyominoes_to_file(input_filename, output_filename=None):
    if output_filename is None:
        output_filename = 'composed-polyomino.txt'
    f = open(output_filename, 'w')
    polyominoes = list(composed_polyominoes_iterator_from_file(input_filename))
    for (word,polyomino) in polyominoes:
        word = ''.join(map(str, word))
        f.write('%s:%s\n'%(word,polyomino))
    f.close()

# ----------- #
# Tikz output #
# ----------- #

def tikz(word, polyomino):
    s = ''
    for c in polyomino:
        cp = (c[0] + 1, c[1] + 1)
        s += '  \\draw[cell] %s rectangle %s;\n' % (c, cp)
    s += '  \\draw[boundary] %s'%str(min(polyomino, key=lambda c: (c[1] != 0, c[0])))
    for l in word:
        if l == 0:
            s += ' -- ++ (1,0)'
        elif l == 1:
            s += ' -- ++ (0,1)'
        elif l == 2:
            s += ' -- ++ (-1,0)'
        elif l == 3:
            s += ' -- ++ (0,-1)'
    s += ';\n'
    return s

def tikz_polyomino_matrix(polyominoes, num_cols=10, col_width=10, row_height=10):
    s = '\\begin{tikzpicture}[cell/.style={fill=black!20, densely dotted}, boundary/.style={very thick}]\n'
    i = 0
    for (boundary_word,polyomino) in polyominoes:
        x = (i % num_cols) * col_width
        y = -((i / num_cols) * row_height)
        w = max(map(lambda c: c[0], polyomino))
        h = max(map(lambda c: c[1], polyomino))
        x += (col_width - w) / 2.0
        y += (row_height - h) / 2.0
        s += '  \\begin{scope}[xshift=%scm,yshift=%scm]\n'%(x,y)
        s += tikz(boundary_word,polyomino)
        s += '  \\end{scope}\n'
        i += 1
    s += '\\end{tikzpicture}\n'
    return s

# ---- #
# Main #
# ---- #

from argparse import ArgumentParser
import os
import subprocess

# Arguments parsing
parser = ArgumentParser()
parser.add_argument('output_filename', help='the name of the output file to which polyominoes are written')
parser.add_argument('--max_cells', type=int, help='the maximum number of cells in the polyomino', default=None)
parser.add_argument('--input_filename', help='the name of the input file from which polyominoes are obtained', default=None)
parser.add_argument('--type', help='the type of polyominoes to be computed (free or composed)', default='free')
parser.add_argument('--tikz', help='write to tex file for tikz output', action='store_true')
parser.add_argument('--num_cols', help='the number of columns in the tikz matrix', type=int, default=10)
parser.add_argument('--col_width', help='the width of each column in the tikz matrix', type=int, default=10)
parser.add_argument('--row_height', help='the height of each row in the tikz matrix', type=int, default=10)
args = parser.parse_args()

if args.tikz:
    polyominoes = polyominoes_iterator_from_file(args.input_filename)
    with open(args.output_filename, 'w') as f:
        f.write(tikz_polyomino_matrix(polyominoes,\
                                      num_cols=args.num_cols,\
                                      col_width=args.col_width,\
                                      row_height=args.row_height))
elif args.type.lower() == 'free':
    if args.max_cells is None:
        print 'You must specify the max number of cells'
    else:
        write_free_polyominoes_to_file(args.max_cells,\
                                       args.output_filename,\
                                       sorted=True)
elif args.type.lower() == 'composed':
    write_composed_polyominoes_to_file(args.input_filename,\
                                       args.output_filename)
