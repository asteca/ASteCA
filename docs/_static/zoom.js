// Click on an image in the page content to enlarge it in place
document.addEventListener("DOMContentLoaded", () => {
  mediumZoom("article img", { background: "rgba(0, 0, 0, 0.8)", margin: 24 });
  // Sphinx wraps scaled images (':scale:', ':width:', ':height:') in a link to the
  // image file: block it, so that clicking zooms the image instead of opening it
  document.querySelectorAll("article a.image-reference").forEach((a) => {
    a.addEventListener("click", (e) => e.preventDefault());
  });
});
