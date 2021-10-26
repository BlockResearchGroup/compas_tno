import compas_rhino
from compas_tno.shapes import Shape
from compas_rhino.geometry import RhinoMesh

answer = compas_rhino.rs.GetString("Select a method to input the shape. From", "Cancel", ["IntraExtrados", "IntraExtradosMiddle", "Middle", "Cancel"])

warning = 'Warning: Select only one mesh.'

middle = None

if answer == "IntraExtrados" or answer == "IntraExtradosMiddle":
    selection = compas_rhino.select_meshes(message='Select Intrados Mesh')
    if not selection:
        print('error')
    if len(selection) > 1:
        compas_rhino.display_message(warning)
    guid_intrados = selection[0]
    intrados_RM = RhinoMesh.from_guid(guid_intrados)
    intrados = intrados_RM.to_compas()
    
    compas_rhino.rs.UnselectObjects(selection)

    selection = compas_rhino.select_meshes(message='Select Extrados Mesh')
    if not selection:
        print('error')
    if len(selection) > 1:
        compas_rhino.display_message(warning)
    guid_extrados = selection[0]
    print(guid_extrados)
    extrados_RM = RhinoMesh.from_guid(guid_extrados)
    extrados = extrados_RM.to_compas()

if answer == "Middle" or answer == "IntraExtradosMiddle":
    guid_middle = compas_rhino.select_mesh(message='Select Intrados Mesh')
    if not guid_middle:
        print('error')
    middle_RM = RhinoMesh.from_guid(guid_middle)
    middle = middle_RM.to_compas()

thk = compas_rhino.rs.GetReal("Input thickness value (estimate)", 0.50)
    
shape = Shape.from_meshes(intrados, extrados, middle, data={'type': 'general', 't': 0.0, 'thk': thk})

shape.to_json('/Users/mricardo/compas_dev/me/plugin/test/test.json')

shape = Shape.from_json('/Users/mricardo/compas_dev/me/plugin/test/test.json')
