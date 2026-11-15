# VS Code setup help

## I cannot find the QPP folder

In the GHCi session from the installation check, enter `:! pwd`. This prints
the full folder path without leaving GHCi. In VS Code, choose **File > Open
Folder**. Paste the path into the dialog's path field. If no path field is
visible, press **Command-Shift-G** on macOS or **Ctrl-L** on Linux to reveal
it. Press Enter to navigate to QPP, then confirm the folder selection.

## QPP settings are not active

Open the cloned QPP folder with **File > Open Folder**. Its file list in
**View > Explorer** should contain `qpp.cabal`, `src` and `.vscode`.
If the window says **Restricted Mode**, click that label and trust this
course folder. On Windows, use a window marked **WSL: Ubuntu**.

## A Unicode shortcut does not expand

Check that the language label at the bottom right says **Haskell**.
The project's `.vscode` folder must contain both `settings.json` and
`qpp.code-snippets`. These files ship with QPP; an older clone may lack them.
In **View > Command Palette**, choose **Preferences: Open Workspace Settings**.
Search `@id:editor.tabCompletion` and check that its value is `onlySnippets`.
The shortcuts work independently of Haskell Language Server.

If input still fails, ask for setup help and report which name you typed,
what appeared after Tab, and whether the checks above passed.

## ⊕ and ⊗ look alike

QPP's settings choose the font JetBrains Mono, which draws `⊕` and `⊗`
larger than VS Code's default font. If the two still look alike, the font
is not installed, and VS Code has quietly used its default font instead.
Install the font:

- macOS: `brew install --cask font-jetbrains-mono`, or download it from
  <https://www.jetbrains.com/lp/mono/> and double-click the font files.
- Ubuntu or Debian: `sudo apt install fonts-jetbrains-mono`.
- Windows: download it from <https://www.jetbrains.com/lp/mono/>, open the
  archive, select the files in `fonts/ttf`, right-click and choose
  **Install**. Install it in Windows, not inside WSL: the VS Code window
  draws its text on the Windows side.

Then quit VS Code and start it again. The font applies only in windows
opened on the QPP folder.

## Check type information

Haskell Language Server supplies types and error messages. It may already
be installed, so the extension need not offer a download. The check below
confirms that type information is available.

In **View > Explorer**, open `src/QPP/Syntax.hs`. Use **Edit > Find** to
locate `maximum $ map`, close Find and hover over `map` on that line.
A popup beginning `map ::` confirms that type information is working;
you do not need to interpret the type for this check.

If nothing appears, allow a few minutes for the first load, move the
pointer away and hover again. If it still fails, open **View > Output**
and select **Haskell (QPP)**. Report the missing popup and the last output
messages when asking for setup help.
