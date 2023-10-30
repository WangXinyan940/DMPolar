import subprocess


def run_command(command: str, **kwargs) -> subprocess.CompletedProcess:
    command_list = command.split()
    return subprocess.run(command_list, **kwargs)