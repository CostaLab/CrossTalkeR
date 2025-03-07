document.addEventListener("DOMContentLoaded", function () {
    let menu = document.querySelector('.navbar-nav .dropdown-menu');
  
    // Create a nested submenu structure
    let submenuItem = document.createElement('li');
    submenuItem.className = "dropdown-submenu";
    submenuItem.innerHTML = `
      <a class="dropdown-item dropdown-toggle" href="#">Example: Bone Marrow Fibrosis in Human (Leimkuhler et. al., 2021)</a>
      <ul class="dropdown-menu">
        <li><a class="dropdown-item" href="articles/HumanFibrosis.html">CrossTalkeR Analysis</a></li>
        <li><a class="dropdown-item" href="articles/LR2TF_analysis.html">Intracellular Communication Analysis</a></li>
        <li><a class="dropdown-item" href="articles/ProgenyLRExample.html">Ligand-Receptor Pathway Enrichment with Progeny</a></li>
      </ul>
    `;
  
    // Find the submenu1 and append the nested menu
    //let submenu1 = menu.querySelector(('h6.dropdown-header[data-toc-skip=Example: Bone Marrow Fibrosis in Human (Leimkuhler et. al., 2021)]')).parentElement;
    // let submenu1 = document.evaluate(
    //     '//h6[contains(text(), "Example: Bone Marrow Fibrosis")]',
    //     document,
    //     null,
    //     XPathResult.FIRST_ORDERED_NODE_TYPE,
    //     null
    //   ).singleNodeValue;
    // submenu1.after(submenuItem);
    let submenu1 = menu.querySelector("a[href='articles/run_liana.html']").parentElement;
    submenu1.appendChild(submenuItem);
    
    // Enable Bootstrap dropdown behavior
    document.querySelectorAll('.dropdown-submenu .dropdown-toggle').forEach((item) => {
      item.addEventListener('click', function (e) {
        e.stopPropagation();
        this.nextElementSibling.classList.toggle('show');
      });
    });
  });