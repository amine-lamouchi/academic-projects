import weakref
import math

GROUP = [
  "amine.lamouchi@polytechnique.edu",
  "giuseppe.cognata@polytechnique.edu",
]

HC = weakref.WeakValueDictionary()

def hc(s):
    return HC.setdefault(s, s)

class Universe:
    def round(self):
        """Compute (in place) the next generation of the universe"""
        raise NotImplementedError

    def get(self, i, j):
        """Returns the state of the cell at coordinates (ij[0], ij[1])"""
        raise NotImplementedError

    def rounds(self, n):
        """Compute (in place) the n-th next generation of the universe"""
        for _i in range(n):
            self.round()

class NaiveUniverse(Universe):
    def __init__(self, n, m, cells):
        self.n = n 
        self.m = m
        self.cells =cells

    def round(self):
        old_cells = self.cells
        n = self.n
        m = self.m

        old_cells = [[None]*m] + old_cells + [[None]*m]
        for i in range(n+2):
                old_cells[i] = [None]+ old_cells[i] +[None]
        
        for i in range(1,n+1):
            for j in range(1,m+1):
                neighbours_list = []
                for k1 in range(-1,2):
                    for k2 in range(-1,2):
                        neighbours_list.append(old_cells[i+k1][j+k2])
                if old_cells[i][j] is False:
                    if neighbours_list.count(True) == 3:
                        self.cells[i-1][j-1] = True
                else:
                    if neighbours_list.count(True) != 3 and neighbours_list.count(True) != 4:
                        self.cells[i-1][j-1] = False
                        
    def get(self, i, j):
        return self.cells[i][j]

class AbstractNode:
    BIG = True

    def __init__(self):
        self._hash  = None
        self._cache = None

    def __hash__(self):
        if self._hash is None:
            self._hash = (
                self.population,
                self.level     ,
                self.nw        ,
                self.ne        ,
                self.sw        ,
                self.se        ,
            )
            self._hash = hash(self._hash)
        return self._hash
        
    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, AbstractNode):
            return False
        return \
            self.level      == other.level      and \
            self.population == other.population and \
            self.nw         is other.nw         and \
            self.ne         is other.ne         and \
            self.sw         is other.sw         and \
            self.se         is other.se

    @property
    def level(self):
        """Level of this node"""
        raise NotImplementedError()

    @property
    def population(self):
        """Total population of the area"""
        raise NotImplementedError()

    nw      = property(lambda self : None)
    ne      = property(lambda self : None)
    sw      = property(lambda self : None)
    se      = property(lambda self : None)

    @staticmethod
    def zero(k):
        if k == 0 :
            res = AbstractNode.cell(False)
        else:
            zero = AbstractNode.zero(k-1)
            res = AbstractNode.node(nw = zero, ne = zero, sw = zero, se = zero)
        return res

    def extend(self):
        k = self.level
        zero = AbstractNode.zero(k)
        zero_internal_node = AbstractNode.zero(k-1)
        if k == 0 :
            return AbstractNode.node(nw = zero, ne = self, sw = zero, se= zero)
        nw1 = AbstractNode.node(nw = zero_internal_node   , ne = zero_internal_node  , sw=zero_internal_node , se = self.nw           )
        ne1 = AbstractNode.node(nw = zero_internal_node   , ne = zero_internal_node  , sw=self.ne            , se = zero_internal_node)
        sw1 = AbstractNode.node(nw = zero_internal_node   , ne = self.sw             , sw=zero_internal_node , se = zero_internal_node) 
        se1 = AbstractNode.node(nw = self.se              , ne = zero_internal_node  , sw=zero_internal_node, se = zero_internal_node )
        res = AbstractNode.node(nw = nw1, ne = ne1, sw = sw1, se = se1)
        return res

    def subuniverse(self):
        Qnw = self.nw
        Qne = self.ne
        Qsw = self.sw
        Qse = self.se
        Qtc = AbstractNode.node(nw = Qnw.ne, ne = Qne.nw, sw = Qnw.se, se = Qne.sw)
        Qbc = AbstractNode.node(nw = Qsw.ne, ne = Qse.nw, sw = Qsw.se, se = Qse.sw)
        Qcc = AbstractNode.node(nw = Qtc.sw, ne = Qtc.se, sw = Qbc.nw, se = Qbc.ne)
        Qcl = AbstractNode.node(nw = Qnw.sw, ne = Qcc.nw, sw = Qsw.nw, se = Qcc.sw)
        Qcr = AbstractNode.node(nw = Qcc.ne, ne = Qne.se, sw = Qcc.se, se = Qse.ne)
        d = {'nw': Qnw, 'ne': Qne, 'sw': Qsw, 'se': Qse, 'tc': Qtc, 'bc': Qbc, 'cc': Qcc, 'cl': Qcl, 'cr': Qcr}
        return d

    def central_part(self):
        Qnw = self.nw
        Qne = self.ne
        Qsw = self.sw
        Qse = self.se
        ne1 = Qne.sw
        nw1 = Qnw.se
        sw1 = Qsw.ne
        se1 = Qse.nw
        return AbstractNode.node(nw = nw1, ne = ne1, sw = sw1, se = se1)


    def forward(self, l = None):
        l = self.level - 2 if l is None else l
        self._cache = {} if self._cache is None else self._cache

        if l in self._cache:
            return self._cache[l]

        if self.population == 0:
            nw1 = self.nw.se
            ne1 = self.ne.sw
            sw1 = self.sw.ne
            se1 = self.se.nw
            c = AbstractNode.node(nw = nw1, ne = ne1, sw = sw1, se = se1)
            self._cache[l] = c
            return c

        if l == self.level - 2:
   
            if self.level < 2:
                return None

            elif self.level == 2:
                cells = [self.nw.nw.alive, self.nw.ne.alive, self.ne.nw.alive, self.ne.ne.alive,
                        self.nw.sw.alive, self.nw.se.alive, self.ne.sw.alive, self.ne.se.alive,
                        self.sw.nw.alive, self.sw.ne.alive, self.se.nw.alive, self.se.ne.alive,
                        self.sw.sw.alive, self.sw.se.alive, self.se.sw.alive, self.se.se.alive]
                w = 0b0
                for k in range(16):
                    if cells[k]:
                        w = w | (1 << (15 - k))

                e_list = 16*[0]
                e_list[15] = 0b0100110000000000
                e_list[12] = 0b0010001100000000
                e_list[0]  = 0b0000000000110010
                e_list[3]  = 0b0000000011000100
                e_list[14] = 0b1010111000000000
                e_list[2]  = 0b0000000011101010
                e_list[8]  = 0b0011001000110000
                e_list[11] = 0b1100010011000000
                e_list[10] = 0b1110101011100000
                e_list[9]  = e_list[10] >> 1
                e_list[6]  = e_list[10] >> 4
                e_list[5]  = e_list[10] >> 5
                e_list[13] = e_list[14] >> 1
                e_list[1]  = e_list[2]  >> 1
                e_list[4]  = e_list[8]  >> 4
                e_list[7]  = e_list[11] >> 4

                alive_neighbours = [0]*16
                for i in range(16):
                    w_prime = w&e_list[i]
                    j = 0
                    while w_prime != 0:
                        j += 1
                        w_prime = w_prime&(w_prime - 1)
                    alive_neighbours[i] = j

                newcells = [False]*16
                for i in range(16):
                    if (not cells[i]) and (alive_neighbours[15 - i] == 3):
                        newcells[i] = True
                    elif cells[i] and (alive_neighbours[15 - i] == 2 or alive_neighbours[15 - i] == 3):
                        newcells[i] = True
                 
                c = AbstractNode.node(nw = AbstractNode.cell(newcells[5]), ne = AbstractNode.cell(newcells[6]), 
                                     sw = AbstractNode.cell(newcells[9]), se = AbstractNode.cell(newcells[10]))
                self._cache[l] = c
                return c

            q = self.subuniverse()
            r = {k : v.forward() for (k,v) in q.items()}
            a = {'nw': AbstractNode.node(nw = r['nw'], ne = r['tc'], sw = r['cl'], se = r['cc']), 
                'ne':  AbstractNode.node(nw = r['tc'], ne = r['ne'], sw = r['cc'], se = r['cr']),
                'sw' : AbstractNode.node(nw = r['cl'], ne = r['cc'], sw = r['sw'], se = r['bc']),
                'se' : AbstractNode.node(nw = r['cc'], ne = r['cr'], sw = r['bc'], se = r['se'])}
            cnw = a['nw'].forward()
            cne = a['ne'].forward()
            csw = a['sw'].forward()
            cse = a['se'].forward()
            c = AbstractNode.node(nw = cnw, ne = cne, sw = csw, se = cse)
            self._cache[l] = c
            return c

        else:
                q = self.subuniverse()
                r = {k : v.central_part() for (k,v) in q.items()}
                a = {'nw': AbstractNode.node(nw = r['nw'], ne = r['tc'], sw = r['cl'], se = r['cc']), 
                    'ne':  AbstractNode.node(nw = r['tc'], ne = r['ne'], sw = r['cc'], se = r['cr']),
                    'sw' : AbstractNode.node(nw = r['cl'], ne = r['cc'], sw = r['sw'], se = r['bc']),
                    'se' : AbstractNode.node(nw = r['cc'], ne = r['cr'], sw = r['bc'], se = r['se'])}
                cnw = a['nw'].forward(l)
                cne = a['ne'].forward(l)
                csw = a['sw'].forward(l)
                cse = a['se'].forward(l)
                c = AbstractNode.node(nw = cnw, ne = cne, sw = csw, se = cse)
                self._cache[l] = c
                return c



    @staticmethod
    def canon(node):
        return hc(node)

    @staticmethod
    def cell(alive):
        return AbstractNode.canon(CellNode(alive))

    @staticmethod
    def node(nw, ne, sw, se):
        return AbstractNode.canon(Node(nw, ne, sw, se))
        
class CellNode(AbstractNode):
    def __init__(self, alive):
        super().__init__()

        self._alive = bool(alive)

    level      = property(lambda self : 0)
    population = property(lambda self : int(self._alive))
    alive      = property(lambda self : self._alive)

class Node(AbstractNode):
    def __init__(self, nw, ne, sw, se):
        super().__init__()

        self._level      = 1 + nw.level
        self._population =  \
            nw.population + \
            ne.population + \
            sw.population + \
            se.population
        self._nw = nw
        self._ne = ne
        self._sw = sw
        self._se = se

    level      = property(lambda self : self._level)
    population = property(lambda self : self._population)

    nw = property(lambda self : self._nw)
    ne = property(lambda self : self._ne)
    sw = property(lambda self : self._sw)
    se = property(lambda self : self._se)

class HashLifeUniverse(Universe):
    def __init__(self, *args):
        if len(args) == 1:
            self._root = args[0]
        else:
            self._root = HashLifeUniverse.load(*args)

        self._generation = 0

    @staticmethod
    def load(n, m, cells):
        level = math.ceil(math.log(max(1, n, m), 2))

        mkcell = getattr(AbstractNode, 'cell', CellNode)
        mknode = getattr(AbstractNode, 'node', Node    )

        def get(i, j):
            i, j = i + n // 2, j + m // 2
            return \
                i in range(n) and \
                j in range(m) and \
                cells[i][j]
                
        def create(i, j, level):
            if level == 0:
                return mkcell(get (i, j))

            noffset = 1 if level < 2 else 1 << (level - 2)
            poffset = 0 if level < 2 else 1 << (level - 2)

            nw = create(i-noffset, j+poffset, level - 1)
            sw = create(i-noffset, j-noffset, level - 1)
            ne = create(i+poffset, j+poffset, level - 1)
            se = create(i+poffset, j-noffset, level - 1)

            return mknode(nw=nw, ne=ne, sw=sw, se=se)
                
        return create(0, 0, level)

    def get(self, i, j):

        k = self._root.level 
        if (i < -2**(k-1)) or (i >= 2**(k-1)) or (j < -2**(k-1)) or (j >= 2**(k-1)):
            return False

        if k == 1:
            if i == 0 and j == 0:
                return self._root.ne.alive
            elif i == -1 and j == 0:
                return self._root.nw.alive
            elif i == -1 and j == -1:
                return self._root.sw.alive
            elif i == 0 and j == -1:
                return self._root.se.alive

        if i >= 0 and j >= 0:
            res = HashLifeUniverse(self._root.ne).get(i - 2**(k-2), j - 2**(k-2))
        elif i < 0 and j >= 0:
            res = HashLifeUniverse(self._root.nw).get(i + 2**(k-2), j - 2**(k-2))
        elif i < 0 and j < 0:
            res = HashLifeUniverse(self._root.sw).get(i + 2**(k-2), j + 2**(k-2))
        elif i >= 0 and j < 0:
            res = HashLifeUniverse(self._root.se).get(i - 2**(k-2), j + 2**(k-2))
        return res

    def peripheral_pop(self):
        Qne = self._root.ne
        Qnw = self._root.nw
        Qsw = self._root.sw
        Qse = self._root.se
        res = Qne.ne.population + Qne.nw.population + Qne.se.population +\
              Qnw.ne.population + Qnw.nw.population + Qnw.sw.population +\
              Qsw.nw.population + Qsw.sw.population + Qsw.se.population +\
              Qse.ne.population + Qne.se.population + Qse.sw.population 
        return res

    def extend(self, k):
        while self._root.level < max(2, k) and self.peripheral_pop != 0:
            self._root = self._root.extend()
        self._root = self._root.extend()

    def rounds(self, n):
        binary = bin(n)[2:]
        for i in range(0, len(binary)):
            if binary[len(binary)-i-1] == '1':
                self.extend(i+2)
                self._root = self._root.forward(i)
        self._generation += n

    def round(self):
        return self.rounds(1)

    @property
    def root(self):
        return self._root
        
    @property
    def generation(self):
        return self._generation