document.addEventListener("DOMContentLoaded", function () {
    console.log("extra.js loaded!"); // Debugging check
    
    let menu = document.querySelector('.navbar-nav .dropdown-menu');
  
    // Create a nested submenu structure
    let submenuItem = document.createElement('li');
    submenuItem.className = "dropdown-submenu";
    submenuItem.innerHTML = `
      <a class="dropdown-item dropdown-toggle" href="#">Example: Bone Marrow Fibrosis in Human (Leimkuhler et. al., 2021)</a>
      <ul class="dropdown-menu">
        <li><a class="dropdown-item subsub" href="//costalab.github.io/CrossTalkeR/articles/HumanFibrosis.html">CrossTalkeR Analysis</a></li>
        <li><a class="dropdown-item subsub" href="//costalab.github.io/CrossTalkeR/articles/LR2TF_analysis.html">Intracellular Communication Analysis</a></li>
        <li><a class="dropdown-item subsub" href="//costalab.github.io/CrossTalkeR/articles/ProgenyLRExample.html">Ligand-Receptor Pathway Enrichment with Progeny</a></li>
      </ul>
    `;

    //let submenu1 = menu.querySelector("a[href='articles/run_liana.html']")?.parentElement;
  
    let submenu1 = Array.from(menu.querySelectorAll("a"))
      .find(link => 
      link.getAttribute("href") === "articles/run_liana.html" || 
      link.getAttribute("href") === "../articles/run_liana.html"
      )?.parentElement;

    if (!submenu1) {
      console.error("Submenu1 not found!");
      return;
    }

    submenu1.appendChild(submenuItem);

    // Enable Bootstrap dropdown behavior
    document.querySelectorAll('.dropdown-submenu .dropdown-toggle').forEach((item) => {
      item.addEventListener('click', function (e) {
        e.stopPropagation();
        this.nextElementSibling.classList.toggle('show');
      });
    });
  });

  