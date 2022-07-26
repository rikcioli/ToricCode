from typing import Dict

# TODO: DOCUMENTATION
class LatticeDrawer(object):

    default_char_table = dict( hlink='--', vlink='|', site='+')

    def __init__(self, lattice, char_table: Dict=None):
        self.lattice = lattice
        self.Lx      = self.lattice.Lx
        self.pbc_x   = self.lattice.pbc_x
        self.char_table = char_table if char_table is not None \
                            else self.default_char_table

    def draw_row_links(self, y: int) -> str:
        """Draw the horizontal links of a given row"""
        hlink = self.char_table['hlink']
        site  = self.char_table['site']
        blank  = ' ' * len(hlink)
        bigblank = blank * 2 + '  '
        bulk, edges = self.lattice.hlinks(y, from_zero=True)
        left  = blank + f'{edges[0]:>2d}' + hlink if self.pbc_x else bigblank
        right = hlink + f'{edges[1]:<2d}' + blank if self.pbc_x else bigblank
        left  = left  + site + hlink
        right = hlink + site + right
        inner = (hlink + site + hlink).join(f'{link:2d}' for link in bulk)
        return left + inner + right + '\n'

    def draw_vert_empty(self) -> str:
        """Draw empty vertical spaces"""
        blank = ' ' * (2 * len(self.char_table['hlink']) + 2)
        vlink = self.char_table['vlink']
        return vlink.join(blank for _ in range(self.Lx+1)) + '\n'

    def draw_vert_links(self, y: int) -> str:
        """Draw the vertical links at a given y"""
        links = self.lattice.vlinks(y, from_zero=True)
        blank = ' ' * (2 * len(self.char_table['hlink']) + 2)
        return blank + blank[1:].join(f'{link:<2d}' for link in links) + '\n'

    def draw_lattice(self) -> str:
        """Draw the lattice"""
        picture = ''
        for y in range(self.lattice.Ly-1, 0, -1):
            picture += self.draw_row_links(y) \
                 + self.draw_vert_empty() \
                 + self.draw_vert_links(y-1) \
                 + self.draw_vert_empty()
        picture += self.draw_row_links(0)

        if self.lattice.pbc_y:
            picture = self.draw_vert_links(self.lattice.Ly-1) \
                + self.draw_vert_empty() \
                + picture \
                + self.draw_vert_empty() \
                + self.draw_vert_links(self.lattice.Ly-1)
        return picture
