---
title: "Install & Guide"
permalink: /categories/
layout: tags
author_profile: true
sort_by: date 
sort_order: descending
---

{% assign filtered_posts = site.posts | where_exp: "post", "post.tags contains 'Install' or post.tags contains 'Guide'" %}
{% assign sorted_posts = filtered_posts | sort: 'date' | reverse %}

<ul>
  {% for post in sorted_posts %}
  <li>
    <a href="{{ post.url }}">{{ post.title }}</a> - {{ post.date | date: "%Y-%m-%d" }}
  </li>
  {% endfor %}
</ul>
