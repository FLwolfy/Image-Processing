import cv2
import numpy as np

def warp_a4_from_image(raw_image_path, output_image_path):
    # Step 1: Read the image
    with open(raw_image_path, "rb") as f:
        raw_data = np.fromfile(f, dtype=np.uint8)  # Read raw file as bytes
    img = cv2.imdecode(raw_data, cv2.IMREAD_COLOR)
    
    # Step 2: Convert to grayscale
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Step 3: Edge detection
    edges = cv2.Canny(gray, 50, 150)

    # Step 4: Find contours
    contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)  # Sort by area

    # Step 5: Find the biggest quadrilateral
    a4_contour = None
    for contour in contours:
        approx = cv2.approxPolyDP(contour, 0.02 * cv2.arcLength(contour, True), True)
        if len(approx) == 4:
            a4_contour = approx
            break

    if a4_contour is None:
        raise ValueError("Could not find A4 paper in the image!")

    # Step 6: Get points and warp
    points = a4_contour.reshape(4, 2)
    # Sort the points to [top-left, top-right, bottom-right, bottom-left]
    rect = np.zeros((4, 2), dtype="float32")
    s = points.sum(axis=1)
    rect[0] = points[np.argmin(s)]  # Top-left
    rect[2] = points[np.argmax(s)]  # Bottom-right
    diff = np.diff(points, axis=1)
    rect[1] = points[np.argmin(diff)]  # Top-right
    rect[3] = points[np.argmax(diff)]  # Bottom-left

    # Dimensions for A4 paper
    width = 210 * 2  # mm to pixels (example scaling)
    height = 297 * 2

    dst = np.array([
        [0, 0],
        [width - 1, 0],
        [width - 1, height - 1],
        [0, height - 1]
    ], dtype="float32")

    # Perspective transformation
    matrix = cv2.getPerspectiveTransform(rect, dst)
    warped = cv2.warpPerspective(img, matrix, (width, height))

    # Save the result
    cv2.imwrite(output_image_path, warped)
    print(f"Warped image saved to {output_image_path}")

# Example usage
raw_image_path = "images/Project/paper3.jpg"
output_image_path = "output/paper3.jpg"
warp_a4_from_image(raw_image_path, output_image_path)
