import click

def main():
    pass

app = click.Command()
app.add_command(main)

if __name__ == "__main__":
    app()
