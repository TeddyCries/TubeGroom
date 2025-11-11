"""
TubeGroom: tools for creating, editing, and managing tube-based geometry in Blender.
"""

import importlib
import sys
import bpy
from bpy.app.handlers import persistent
from . import core

bl_info = {
    "name": "TubeGroom",
    "author": "TeddyCry",
    "version": (1, 0, 0),
    "blender": (4, 5, 0),
    "location": "View3D > Sidebar > TubeGroom",
    "description": "Interactive mesh grooming tools",
    "category": "Mesh",
}

MODULES = ["prop", "ui"]

_loaded = []

@persistent
def _clear_state(_):
    if hasattr(core, "reset_all_data"):
        try:
            core.reset_all_data(clear_history=True)
        except Exception:
            pass

def register():
    global _loaded
    _loaded = []
    for name in MODULES:
        full = f"{__package__}.{name}"
        mod = importlib.reload(sys.modules[full]) if full in sys.modules else importlib.import_module(f".{name}", __package__)
        _loaded.append(mod)
        if hasattr(mod, "register"):
            mod.register()
    if _clear_state not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(_clear_state)

def unregister():
    try:
        bpy.app.handlers.load_post.remove(_clear_state)
    except Exception:
        pass
    for mod in reversed(_loaded):
        if hasattr(mod, "unregister"):
            try:
                mod.unregister()
            except Exception:
                pass
