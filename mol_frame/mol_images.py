import base64
from io import BytesIO as IO

from PIL import Image, ImageChops

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18


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


def b64_mol(mol, size=300):
    img_file = IO()
    if isinstance(mol, list):
        img = autocrop(Draw.MolsToGridImage(mol, size=(size, size)))
    else:
        try:
            img = autocrop(Draw.MolToImage(mol, size=(size, size)))
        except UnicodeEncodeError:
            print(Chem.MolToSmiles(mol))
            mol = Chem.MolFromSmiles("*")
            img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img = make_transparent(img)
    img.save(img_file, format='PNG')
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


def mol_img_tag(mol, options=None):
    tag = """<img {} src="data:image/png;base64,{}" alt="Mol"/>"""
    if options is None:
        options = ""
    img_tag = tag.format(options, b64_mol(mol))
    return img_tag
