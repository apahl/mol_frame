import base64
from itertools import chain
from io import BytesIO as IO

import numpy as np

from PIL import Image, ImageChops

from rdkit.Chem import AllChem as Chem
from rdkit.Chem.rdCoordGen import AddCoords  # New coord. generation
from rdkit.Chem import Draw

from mol_frame import tools as mft

try:
    Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
    Draw.DrawingOptions.atomLabelFontSize = 18
except KeyError:  # Font "DejaVu Sans" is not available
    pass


config = mft.load_config()

USE_RDKIT_NEW_COORD = config["Options"].get("UseNewRdkitCoord", True)

# try:
#     # Try to import Avalon so it can be used for generation of 2d coordinates.
#     from rdkit.Avalon import pyAvalonTools as pyAv
#     USE_AVALON_2D = True
# except ImportError:
#     print("* Avalon not available. Using RDKit for 2d coordinate generation.")
#     USE_AVALON_2D = False


def rescale(mol, f=1.4):
    tm = np.zeros((4, 4), np.double)
    for i in range(3):
        tm[i, i] = f
    tm[3, 3] = 1.0
    Chem.TransformMol(mol, tm)


def add_coords(mol, force=False):
    """Check if a mol has 2D coordinates and if not, calculate them."""
    if not force:
        try:
            mol.GetConformer()
        except ValueError:
            force = True  # no 2D coords... calculate them

    if force:
        if USE_RDKIT_NEW_COORD and mol.GetNumAtoms() <= 75:
            AddCoords(mol)
            rescale(mol, f=1.4)
        else:
            mol.Compute2DCoords()


def make_transparent(img):
    img = img.convert("RGBA")
    pixdata = img.load()
    width, height = img.size
    for y in range(height):
        for x in range(width):
            if pixdata[x, y] == (255, 255, 255, 255):
                pixdata[x, y] = (255, 255, 255, 0)
    return img


def autocrop(im, bgcolor="white"):
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = Image.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    return None  # no contents


def mol_img_file(mol, size=300, hlsss=None):
    img_file = IO()
    if isinstance(mol, list):
        [add_coords(x) for x in mol]
        img = autocrop(Draw.MolsToGridImage(mol, size=(size, size)))
    else:
        if isinstance(mol, str):  # convert from Smiles on-the-fly, when necessary
            mol = Chem.MolFromSmiles(mol)
        if mol is None:
            mol = Chem.MolFromSmiles("*")
        if hlsss is not None:
            if isinstance(hlsss, str):
                hlsss = hlsss.split(",")
                atoms = set()
                for smi in hlsss:
                    m = Chem.MolFromSmiles(smi)
                    if m:
                        matches = list(chain(*mol.GetSubstructMatches(m)))
                    else:
                        matches = []
                    if len(matches) > 0:
                        atoms = atoms.union(set(matches))
            atoms = list(atoms)
        else:
            atoms = []
        try:
            add_coords(mol)
            img = autocrop(
                Draw.MolToImage(mol, size=(size, size), highlightAtoms=atoms)
            )
        except UnicodeEncodeError:
            print(Chem.MolToSmiles(mol))
            mol = Chem.MolFromSmiles("*")
            img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img = make_transparent(img)
    img.save(img_file, format="PNG")
    val = img_file.getvalue()
    img_file.close()
    return val


def b64_mol(mol, size=300, hlsss=None):
    img_file = mol_img_file(mol, size=size, hlsss=hlsss)
    b64 = base64.b64encode(img_file)
    b64 = b64.decode()
    return b64


def mol_img_tag(mol, size=300, options=None):
    tag = """<img {} src="data:image/png;base64,{}" alt="Mol"/>"""
    if options is None:
        options = ""
    img_tag = tag.format(options, b64_mol(mol, size=size))
    return img_tag
