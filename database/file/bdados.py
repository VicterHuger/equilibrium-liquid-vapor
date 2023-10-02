def data_base(system_id: str):
    """A function to read database as file"""

    with open("./DADOS.DAT", encoding="utf-8") as bdados_file:
        file_content = bdados_file.read()
        print(file_content)
