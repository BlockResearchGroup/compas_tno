from compas.geometry import Point, Polyline, Bezier
from compas.colors import Color
from compas_view2.app import App

# future

# curve = Bezier([[0, 0, 0], [3, 6, 0], [5, -3, 0], [10, 0, 0]])

# viewer = App(viewmode="shaded", enable_sidebar=True, width=1600, height=900)
# viewer.view.camera.tx = -5.0
# viewer.view.camera.rz = 0
# viewer.view.camera.rx = -20

# pointobj = viewer.add(Point(* curve.point(0)), size=20, color=(1, 0, 0))
# curveobj = viewer.add(Polyline(curve.locus()), linewidth=2)


# @viewer.checkbox(text="Show Point", checked=True)
# def check(checked):
#     pointobj.is_visible = checked
#     viewer.view.update()


# @viewer.slider(title="Slide Point", maxval=100, step=1, bgcolor=Color.white())
# def slide(value):
#     value = value / 100
#     pointobj._data = curve.point(value)
#     pointobj.update()
#     viewer.view.update()


# @viewer.button(text="Reset")
# def click():
#     if viewer.confirm('This will reset the point to parameter t=0.'):
#         pointobj._data = curve.point(0)
#         pointobj.update()
#         slide.value = 0
#         viewer.view.update()


# @viewer.radio(title='Display', items=[
#     {'text': 'Ghosted', 'value': 'ghosted', 'checked': viewer.view.mode == 'ghosted'},
#     {'text': 'Shaded', 'value': 'shaded', 'checked': viewer.view.mode == 'shaded'},
#     {'text': 'Lighted', 'value': 'lighted', 'checked': viewer.view.mode == 'lighted'},
#     {'text': 'Wireframe', 'value': 'wireframe', 'checked': viewer.view.mode == 'wireframe'}
# ])
# def select1(value):
#     viewer.view.mode = value
#     viewer.view.update()


# @viewer.select(items=[
#     {'text': 'Item 1'},
#     {'text': 'Item 2'},
#     {'text': 'Item 3'},
#     {'text': 'Item 4'}
# ])
# def select2(index, text):
#     viewer.info(f"You selected item '{index}' with text '{text}'")


# viewer.run()
