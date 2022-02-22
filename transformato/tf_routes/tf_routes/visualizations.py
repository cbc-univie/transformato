from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.ipython_useSVG = False
import logging
from colorsys import hls_to_rgb
from copy import deepcopy

from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D

logger = logging.getLogger(__name__)


def _show_common_core_gradient(
    mol, highlight, mutationl, percomponent=False, numbers=False
):

    """
    Draw the molecule with color gradient for mutation order

    highlight: list indicating the ccore, fct highlights the common core in black
    mutationl: list indicating the mutation order, fct draws color gradient

    percompoentn: boolean; one color gradient for each component (False) or one for all dummy atoms
    - currently, the components are processed successively, so perhaps True is better suited to illustrate the mutation order

    numbers: boolean; show also numbers showing the mutation order
    """

    # https://rdkit.blogspot.com/2015/02/new-drawing-code.html

    mol = deepcopy(mol)
    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(800, 800)

    # new
    drawer.drawOptions().setHighlightColour((0.8, 0.8, 0.8))

    drawer.SetFontSize(0.3)

    opts = drawer.drawOptions()

    for i in mol.GetAtoms():
        opts.atomLabels[i.GetIdx()] = (
            str(i.GetProp("atom_index")) + ":" + i.GetProp("atom_type")
        )

    colours = [(0.8, 0.0, 0.8), (0, 0.15, 0), (0, 0.8, 0.8), (0, 0, 0.8)]
    atom_cols = {}
    for i, at in enumerate(highlight):
        atom_cols[at] = colours[1]

    if percomponent == False:
        mutationlflat = [item for sublist in mutationl for item in sublist]

        # colour gradient
        rcolours = rainbow_color_gradient(n=len(mutationlflat))

        colours = [(0.8, 0.0, 0.8), (0.8, 0.8, 0), (0, 0.8, 0.8), (0, 0, 0.8)]
        for i, at in enumerate(mutationlflat):
            atom_cols[at] = rcolours[i]

        if numbers == True:
            atomlist = mol.GetAtoms()
            for i, at in enumerate(mutationlflat):
                mnr = i + 1
                lbl = "%.0f" % mnr
                atomlist[at].SetProp("atomNote", lbl)

        mutationlflat = mutationlflat + list(highlight)

        drawer.DrawMolecule(
            mol, highlightAtoms=mutationlflat, highlightAtomColors=atom_cols
        )
        Draw.DrawingOptions.includeAtomNumbers = False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")

        img = Draw.MolToImage(
            mol,
            highlightAtoms=mutationlflat,
            highlightAtomColors=atom_cols,
            useSVG=False,
        )

    else:
        mutationlflat = [item for sublist in mutationl for item in sublist]

        colours = [(0.8, 0.0, 0.8), (0.8, 0.8, 0), (0, 0.8, 0.8), (0, 0, 0.8)]

        componentcount = 0
        for mutcomp in mutationl:

            rcolours = rainbow_color_gradient(n=len(mutcomp))
            for i, at in enumerate(mutcomp):

                atom_cols[at] = rcolours[i]

            componentcount = componentcount + 1

        if numbers == True:
            atomlist = mol.GetAtoms()
            for i, at in enumerate(mutationlflat):
                mnr = i + 1
                lbl = "%.0f" % mnr
                atomlist[at].SetProp("atomNote", lbl)

        mutationlflat = mutationlflat + list(highlight)

        drawer.DrawMolecule(
            mol, highlightAtoms=mutationlflat, highlightAtomColors=atom_cols
        )
        Draw.DrawingOptions.includeAtomNumbers = False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")

        img = Draw.MolToImage(
            mol,
            highlightAtoms=mutationlflat,
            highlightAtomColors=atom_cols,
            useSVG=False,
        )

    return svg


def _show_common_core_gradient_write(mol, highlight, mutationl, percomponent=False):

    """
    almost identical to _show_common_core_gradient_write, used for creating visualizations to write in files

    Draw the molecule with color gradient for mutation order

    highlight: list indicating the ccore, fct highlights the common core in black
    mutationl: list indicating the mutation order, fct draws color gradient

    percompoentn: boolean; one color gradient for each component (False) or one for all dummy atoms
    - currently, the components are processed successively, so perhaps True is better suited to illustrate the mutation order

    numbers: boolean; show also numbers showing the mutation order
    legendtext: string; legend to be placed below molecule
    """

    # https://rdkit.blogspot.com/2015/02/new-drawing-code.html

    mol = deepcopy(mol)
    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)

    drawer.drawOptions().setHighlightColour((0.8, 0.8, 0.8))

    drawer.SetFontSize(0.3)

    opts = drawer.drawOptions()

    for i in mol.GetAtoms():
        opts.atomLabels[i.GetIdx()] = (
            str(i.GetProp("atom_index")) + ":" + i.GetProp("atom_type")
        )

    colours = [(0.8, 0.0, 0.8), (0, 0.15, 0), (0, 0.8, 0.8), (0, 0, 0.8)]
    atom_cols = {}
    for i, at in enumerate(highlight):
        atom_cols[at] = colours[1]

    if percomponent == False:
        mutationlflat = [item for sublist in mutationl for item in sublist]

        # rainbow colors
        rcolours = rainbow_color_gradient(n=len(mutationlflat))

        colours = [(0.8, 0.0, 0.8), (0.8, 0.8, 0), (0, 0.8, 0.8), (0, 0, 0.8)]
        for i, at in enumerate(mutationlflat):
            atom_cols[at] = rcolours[i]

        mutationlflat = mutationlflat + list(highlight)

        drawer.DrawMolecule(
            mol, highlightAtoms=mutationlflat, highlightAtomColors=atom_cols
        )
        Draw.DrawingOptions.includeAtomNumbers = False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")

        img = Draw.MolToImage(
            mol,
            highlightAtoms=mutationlflat,
            highlightAtomColors=atom_cols,
            useSVG=False,
        )

    else:
        mutationlflat = [item for sublist in mutationl for item in sublist]

        colours = [(0.8, 0.0, 0.8), (0.8, 0.8, 0), (0, 0.8, 0.8), (0, 0, 0.8)]

        componentcount = 0
        for mutcomp in mutationl:

            rcolours = rainbow_color_gradient(n=len(mutcomp))
            for i, at in enumerate(mutcomp):

                atom_cols[at] = rcolours[i]

            componentcount = componentcount + 1

        mutationlflat = mutationlflat + list(highlight)

        drawer.DrawMolecule(
            mol, highlightAtoms=mutationlflat, highlightAtomColors=atom_cols
        )
        Draw.DrawingOptions.includeAtomNumbers = False
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")

        img = Draw.MolToImage(
            mol,
            highlightAtoms=mutationlflat,
            highlightAtomColors=atom_cols,
            useSVG=False,
        )

    return svg, img


def animated_visualization_3d_v1(mol, mutationl, ccoremol, hits):
    """
    animated visualization using py3Dmol
    """

    import time

    import py3Dmol

    ccoremol
    ccoretemplate = ccoremol
    AllChem.EmbedMolecule(ccoretemplate)

    m_match = hits

    molsreduced = []
    flatorder = [item for items in mutationl for item in items]

    for i in range(len(flatorder)):

        mol1newembed = AllChem.ConstrainedEmbed(mol, ccoretemplate)

        mol1ed = Chem.EditableMol(mol1newembed)
        currentorderlist = flatorder[0 : i + 1]
        currentorderlist.sort(reverse=True)

        for j in currentorderlist:
            mol1ed.RemoveAtom(j)
        molback = mol1ed.GetMol()

        molsreduced.append(molback)

    ccoreblock = Chem.MolToMolBlock(ccoretemplate, kekulize=False)

    molblocksreduced = []
    for i in molsreduced:
        mblock = Chem.MolToMolBlock(i, kekulize=False)

        molblocksreduced.append(mblock)

    view = py3Dmol.view()

    counter = 0
    for i in molblocksreduced:

        view.removeAllModels()

        view.addModel(i, "mol")
        view.setStyle({"model": 0}, {"stick": {}})

        view.setBackgroundColor("0xeeeeee")

        if counter == 0:
            view.zoomTo()
            counter = counter + 1

        view.update()
        time.sleep(1)


def animated_visualization_3d_v2(mol, mutationl, ccoremol, hits):
    """
    animated visualization using py3Dmol
    common core molecules are visualized as spheres
    """

    import time

    import py3Dmol
    from rdkit.Chem import TemplateAlign

    m_match = hits

    ccoremol
    ccoretemplate = ccoremol
    AllChem.EmbedMolecule(ccoretemplate)

    mol1newembed = AllChem.ConstrainedEmbed(mol, ccoretemplate)

    molsreduced = []

    flatorder = [item for items in mutationl for item in items]

    newhits = mol.GetSubstructMatch(ccoremol)
    mol.GetNumConformers()
    counter = 0
    for i in range(mol.GetNumAtoms()):

        for j in newhits:

            if j == i:

                newpos = ccoremol.GetConformer().GetAtomPosition(counter)
                mol.GetConformer().SetAtomPosition(i, newpos)

                counter = counter + 1

    for i in range(len(flatorder)):

        mol1newembed = mol

        mol1newembed = AllChem.ConstrainedEmbed(mol1newembed, ccoretemplate)

        mol1newembed.GetNumConformers()

        mol1ed = Chem.EditableMol(mol1newembed)
        currentorderlist = flatorder[0 : i + 1]
        currentorderlist.sort(reverse=True)
        for j in currentorderlist:
            mol1ed.RemoveAtom(j)
        molback = mol1ed.GetMol()

        molsreduced.append(molback)

    ccoreblock = Chem.MolToMolBlock(ccoretemplate, kekulize=False)

    molblocksreduced = []
    for i in molsreduced:
        mblock = Chem.MolToMolBlock(i, kekulize=False)

        molblocksreduced.append(mblock)

    view = py3Dmol.view()

    counter = 0
    for i in molblocksreduced:

        view.removeAllModels()

        view.addModel(i, "mol")
        view.setStyle({"model": 0}, {"stick": {}})

        view.addModel(ccoreblock, "mol")
        # view.setStyle({'model':1},{'sphere':{'color':'black', 'radius':1.8}})
        view.setStyle({"model": 1}, {"sphere": {"radius": 0.8}})
        view.setBackgroundColor("0xeeeeee")

        if counter == 0:
            view.zoomTo()
            counter = counter + 1

        view.update()

        time.sleep(1)


def convertgraphics(img1, img2):
    """
    converting files to appropriate format for writing in file
    """

    from io import BytesIO

    from cairosvg import svg2png

    imgn1 = svg2png(img1)
    imgn2 = svg2png(img2)

    import PIL.Image

    imo1 = PIL.Image.open(BytesIO(imgn1))
    imo2 = PIL.Image.open(BytesIO(imgn2))

    return imo1, imo2


def rainbow_color_gradient(n=10, end=1 / 3):
    #  end=n/40
    #  end = 1/3
    if n > 1:
        return [hls_to_rgb(end * i / (n - 1), 0.5, 1) for i in range(n)]
    else:
        return [hls_to_rgb(end * i / (n), 0.5, 1) for i in range(n)]
