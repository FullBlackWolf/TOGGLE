---
title: "Install & Guide"
permalink: /categories/
layout: tags
author_profile: true
sort_by: date 
sort_order: descending
---

<ul>
  {% for post in site.posts %}
    {% if post.tags contains 'Install' or post.tags contains 'Guide' %}
      <li>
        <a href="{{ post.url }}">{{ post.title }}</a> - {{ post.date | date: "%Y-%m-%d" }}
      </li>
    {% endif %}
  {% endfor %}
</ul>
