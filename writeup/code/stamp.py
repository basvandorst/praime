"""Set the red 'last updated' line at the top of README.md to the current Amsterdam time.
Run after every new answer or proof, before the commit."""
import re, subprocess
stamp = subprocess.run(["date", "+%Y-%m-%d %H:%M (%Z)"], env={"TZ": "Europe/Amsterdam"},
                       capture_output=True, text=True, check=True).stdout.strip()
line = f"$\\color{{red}}{{\\textsf{{last updated: {stamp}}}}}$\n"
text = open("README.md", encoding="utf-8").read()
text = re.sub(r"^\$\\color\{red\}.*\n\n", "", text)
open("README.md", "w", encoding="utf-8").write(line + "\n" + text)
print(line.strip())
