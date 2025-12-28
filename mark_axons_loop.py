import os
import sys
import subprocess
from typing import List


# Folder containing your processed *_green_only.tif images
input_dir = r'E:\mGreenLantern_40x\adult_filtered'


def write_macro_for_image(image_path: str, macro_path: str) -> None:
    """Create a Fiji macro that prepares the image for axon marking.

    You draw and save ROI sets / values yourself; this script only opens images.
    """
    image_path_macro = image_path.replace("\\", "/")
    full_base = os.path.splitext(os.path.basename(image_path))[0]
    if full_base.endswith("_green_only"):
        base = full_base[: -len("_green_only")]
    else:
        base = full_base

    macro_content = f"""// Axon marking workflow macro
// Opens one image, sets scale, prepares tools, and shows instructions.

open("{image_path_macro}");

// Optional: set a fixed scale (adjust if needed)
run("Set Scale...", "distance=1 known=0.1625 pixel=1 unit=micron");

// Open Brightness/Contrast so you can adjust visibility if needed
run("Brightness/Contrast...");

// Tool setup: select Straight Line tool and set its line width to 50,
// then switch to Freehand line (under Straight line)
setTool("line");
run("Line Width...", "line=50");
setTool("freehand");

// Open ROI Manager so 'T' adds ROIs there
run("ROI Manager...");
roiManager("Show All");

// Instructions popup: you save files and choose names yourself
waitForUser("Axon marking instructions",
    "IMAGE: {base}\\n\\n" +
    "1. Draw 5-10 axons per image with the Freehand Line tool (under Straight line).\\n" +
    "2. Line Width is already set to 50 via the Straight line tool.\\n" +
    "3. Press 'T' after each axon to add it to ROI Manager.\\n" +
    "4. When you are DONE with this image, manually:\\n" +
    "   - In ROI Manager: 'More >>' -> 'Save...' and choose your ROISet filename.\\n" +
    "   - In Results table: 'File -> Save As...' and choose your values filename.\\n" +
    "\\nAfter saving ROI set and values, close this dialog and then close Fiji.\\n" +
    "The Python script will then open the next image."
);
"""

    with open(macro_path, "w") as f:
        f.write(macro_content)


def launch_fiji_with_macro(macro_path: str) -> subprocess.Popen:
    """Launch Fiji with the given macro and return the process handle."""
    fiji_dir = r"C:\Program Files\fiji-win64\Fiji.app"
    candidates = [
        os.path.join(fiji_dir, "ImageJ-win64.exe"),
        os.path.join(fiji_dir, "Fiji.exe"),
        os.path.join(fiji_dir, "ImageJ.exe"),
    ]

    fiji_path = next((p for p in candidates if os.path.exists(p)), None)
    if not fiji_path:
        raise RuntimeError(f"Fiji executable not found in: {fiji_dir}")

    creation_flags = subprocess.CREATE_NEW_CONSOLE if sys.platform == "win32" else 0
    proc = subprocess.Popen(
        [fiji_path, "-macro", macro_path],
        creationflags=creation_flags,
    )
    return proc


def main() -> None:
    # Find all *_green_only.tif images in the folder
    green_only_files: List[str] = []
    for name in sorted(os.listdir(input_dir)):
        if not name.lower().endswith("_green_only.tif"):
            continue
        green_only_files.append(os.path.join(input_dir, name))

    if not green_only_files:
        print(f"No *_green_only.tif images found in {input_dir}")
        return

    print("Found the following *_green_only.tif images to process:")
    for path in green_only_files:
        print(f" - {path}")

    macro_file = os.path.join(input_dir, "mark_axons_workflow.ijm")

    for idx, img_path in enumerate(green_only_files, start=1):
        print(f"\n=== Image {idx} of {len(green_only_files)} ===")
        print(f"Preparing macro and launching Fiji for: {img_path}")

        # Write macro for this specific image
        write_macro_for_image(img_path, macro_file)

        # Launch Fiji
        try:
            proc = launch_fiji_with_macro(macro_file)
            print("Fiji launched. Follow the on-screen instructions, save ROI set and values manually, then close Fiji.")
        except Exception as e:
            print(f"Could not launch Fiji: {e}")
            print("Please open the image manually in Fiji and follow your workflow.")
            continue

        # Wait for Fiji process to exit (user closes it) before moving to the next image
        try:
            proc.wait()
        except Exception:
            pass

        print(f"Finished with image: {img_path}")

    print("\nFinished. Each *_green_only image has been opened once in Fiji. Any ROI sets / values you saved there are now on disk.")


if __name__ == "__main__":
    main()


