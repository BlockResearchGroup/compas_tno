if __name__ == '__main__':
    import compas_rhino
    from compas_rhino.uninstall import uninstall
    from compas_rhino.uninstall_plugin import uninstall_plugin
    import os

    print("\n", "-"*10, "Removing existing plugins", "-"*10)
    python_plugins_path = compas_rhino._get_python_plugins_path("6.0")
    print("Plugin location: ", python_plugins_path)
    plugins = os.listdir(python_plugins_path)
    for p in plugins:
        uninstall_plugin(p, version="6.0")

    print("\n", "-"*10, "Removing existing packages", "-"*10)
    uninstall()

    print("\n", "-"*10, "Uninstallation is successful", "-"*10)
