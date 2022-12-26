import subprocess

if __name__ == '__main__':
    resolutions = [16, 32, 48, 64, 96, 128, 192, 256]
    for res in resolutions:
        print(res)
        print(subprocess.run(['../cmake-build-release/cages', f'disc_{res}', 'disc_8']).stdout)